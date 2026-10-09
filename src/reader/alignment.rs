//! Indexed BAM/CRAM files: opening, header, index and region queries.
//!
//! A file is opened once ([`AlignmentFile::open`]) to read its header and
//! locate its index; each thread then opens its own [`RegionReader`], so no
//! mutable state is shared.

use std::fs::File;
use std::io::{self, Read};
use std::path::{Path, PathBuf};

use noodles::core::{Position, Region};
use noodles::fasta;
use noodles::sam::{self, alignment::Record};
use noodles_util::alignment::io::IndexedReader;
use noodles_util::alignment::io::indexed_reader::Builder;
use thiserror::Error;

use crate::regions::targets::{Contig, ScanRegion};

/// Alignment file format.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    /// BAM.
    Bam,
    /// CRAM.
    Cram,
}

impl Format {
    /// Detect the format from the first bytes of a file.
    fn detect(path: &Path) -> Result<Self, AlignmentError> {
        let mut magic = [0; 4];
        File::open(path)
            .and_then(|mut file| file.read_exact(&mut magic))
            .map_err(|source| AlignmentError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        match magic {
            [b'C', b'R', b'A', b'M'] => Ok(Self::Cram),
            [0x1f, 0x8b, ..] => Ok(Self::Bam),
            _ => Err(AlignmentError::UnknownFormat {
                path: path.to_path_buf(),
            }),
        }
    }

    /// Index paths to look for, in order.
    fn index_candidates(self, path: &Path) -> Vec<PathBuf> {
        let append = |ext: &str| {
            let mut s = path.as_os_str().to_owned();
            s.push(ext);
            PathBuf::from(s)
        };
        match self {
            Self::Bam => vec![append(".bai"), append(".csi"), path.with_extension("bai")],
            Self::Cram => vec![append(".crai"), path.with_extension("crai")],
        }
    }
}

impl std::fmt::Display for Format {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Bam => "BAM",
            Self::Cram => "CRAM",
        })
    }
}

/// Errors raised while opening or reading an alignment file.
#[derive(Debug, Error)]
pub enum AlignmentError {
    /// The file could not be read.
    #[error("cannot read {path}: {source}")]
    Io {
        /// Path of the file.
        path: PathBuf,
        /// Underlying I/O error.
        source: io::Error,
    },
    /// The file is neither BAM nor CRAM.
    #[error("{path} is neither a BAM nor a CRAM file")]
    UnknownFormat {
        /// Path of the file.
        path: PathBuf,
    },
    /// No index was found next to the file.
    #[error("no index found for {path} (looked for {}); run `samtools index {path}`",
        .candidates.iter().map(|p| p.display().to_string()).collect::<Vec<_>>().join(", "))]
    MissingIndex {
        /// Path of the alignment file.
        path: PathBuf,
        /// Index paths looked for.
        candidates: Vec<PathBuf>,
    },
    /// A CRAM file was given without a reference.
    #[error("{path} is a CRAM file: --reference is required to decode it")]
    MissingReference {
        /// Path of the CRAM file.
        path: PathBuf,
    },
    /// The reference FASTA could not be opened.
    #[error(
        "cannot open the reference {path} (indexed FASTA, run `samtools faidx` if the .fai is missing): {source}"
    )]
    Reference {
        /// Path of the FASTA file.
        path: PathBuf,
        /// Underlying I/O error.
        source: io::Error,
    },
}

/// An indexed BAM or CRAM file.
pub struct AlignmentFile {
    path: PathBuf,
    format: Format,
    index_path: PathBuf,
    header: sam::Header,
    reference: Option<fasta::Repository>,
}

impl std::fmt::Debug for AlignmentFile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AlignmentFile")
            .field("path", &self.path)
            .field("format", &self.format)
            .field("index_path", &self.index_path)
            .finish_non_exhaustive()
    }
}

impl AlignmentFile {
    /// Open a BAM or CRAM file, read its header and locate its index.
    ///
    /// `reference` is the indexed FASTA used to decode a CRAM file; it is
    /// required for CRAM and unused for BAM.
    ///
    /// # Errors
    ///
    /// Fails if the file cannot be read, is neither BAM nor CRAM, has no
    /// index, or is a CRAM file without a usable reference.
    pub fn open(path: &Path, reference: Option<&Path>) -> Result<Self, AlignmentError> {
        let format = Format::detect(path)?;
        let candidates = format.index_candidates(path);
        let index_path = candidates
            .iter()
            .find(|candidate| candidate.is_file())
            .cloned()
            .ok_or_else(|| AlignmentError::MissingIndex {
                path: path.to_path_buf(),
                candidates: candidates.clone(),
            })?;

        let reference = match (format, reference) {
            (Format::Cram, None) => {
                return Err(AlignmentError::MissingReference {
                    path: path.to_path_buf(),
                });
            }
            (Format::Cram, Some(fasta_path)) => Some(open_reference(fasta_path)?),
            (Format::Bam, _) => None,
        };

        let mut file = Self {
            path: path.to_path_buf(),
            format,
            index_path,
            header: sam::Header::default(),
            reference,
        };
        file.header = file
            .build_reader()
            .and_then(|mut reader| reader.read_header())
            .map_err(|source| file.io_error(source))?;
        Ok(file)
    }

    /// Path of the alignment file.
    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// File format.
    #[must_use]
    pub fn format(&self) -> Format {
        self.format
    }

    /// Path of the index in use.
    #[must_use]
    pub fn index_path(&self) -> &Path {
        &self.index_path
    }

    /// SAM header.
    #[must_use]
    pub fn header(&self) -> &sam::Header {
        &self.header
    }

    /// Contigs of the header, in header order.
    #[must_use]
    pub fn contigs(&self) -> Vec<Contig> {
        self.header
            .reference_sequences()
            .iter()
            .map(|(name, contig)| Contig {
                name: name.to_string(),
                length: contig.length().get() as u64,
            })
            .collect()
    }

    /// Sample names (`SM`) of the read groups, deduplicated, in header order.
    #[must_use]
    pub fn samples(&self) -> Vec<String> {
        let mut samples: Vec<String> = Vec::new();
        for read_group in self.header.read_groups().values() {
            let sample = read_group
                .other_fields()
                .get(&sam::header::record::value::map::read_group::tag::SAMPLE);
            if let Some(sample) = sample.map(ToString::to_string)
                && !samples.contains(&sample)
            {
                samples.push(sample);
            }
        }
        samples
    }

    /// Open a reader for region queries; open one per thread.
    ///
    /// # Errors
    ///
    /// Fails if the file or its index cannot be read.
    pub fn reader(&self) -> Result<RegionReader<'_>, AlignmentError> {
        let mut inner = self.build_reader().map_err(|e| self.io_error(e))?;
        // Skip the header to reach the records.
        inner.read_header().map_err(|e| self.io_error(e))?;
        Ok(RegionReader { file: self, inner })
    }

    fn build_reader(&self) -> io::Result<IndexedReader<File>> {
        let builder = Builder::default();
        let mut builder = match self.format {
            Format::Bam if self.index_path.extension().is_some_and(|ext| ext == "csi") => {
                builder.set_index(noodles::csi::fs::read(&self.index_path)?)
            }
            Format::Bam => builder.set_index(noodles::bam::bai::fs::read(&self.index_path)?),
            Format::Cram => builder.set_index(noodles::cram::crai::fs::read(&self.index_path)?),
        };
        if let Some(reference) = &self.reference {
            builder = builder.set_reference_sequence_repository(reference.clone());
        }
        builder.build_from_path(&self.path)
    }

    fn io_error(&self, source: io::Error) -> AlignmentError {
        AlignmentError::Io {
            path: self.path.clone(),
            source,
        }
    }
}

/// Open an indexed FASTA file as a reference sequence repository.
fn open_reference(path: &Path) -> Result<fasta::Repository, AlignmentError> {
    let reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .map_err(|source| AlignmentError::Reference {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(fasta::Repository::new(
        fasta::repository::adapters::IndexedReader::new(reader),
    ))
}

/// A reader of the records overlapping a region.
pub struct RegionReader<'a> {
    file: &'a AlignmentFile,
    inner: IndexedReader<File>,
}

impl RegionReader<'_> {
    /// Records overlapping `region`, in coordinate order. Unmapped reads
    /// placed at the position of their mate are included.
    ///
    /// # Errors
    ///
    /// Fails if the region cannot be queried (unknown contig, unreadable
    /// index); errors while reading records are returned by the iterator.
    pub fn query(
        &mut self,
        region: &ScanRegion,
    ) -> io::Result<impl Iterator<Item = io::Result<Box<dyn Record>>> + '_> {
        let invalid = |e| io::Error::new(io::ErrorKind::InvalidInput, e);
        let start = usize::try_from(region.start + 1)
            .map_err(invalid)
            .and_then(|s| Position::try_from(s).map_err(invalid))?;
        let end = usize::try_from(region.end)
            .map_err(invalid)
            .and_then(|e| Position::try_from(e).map_err(invalid))?;
        let region = Region::new(region.contig.as_str(), start..=end);
        self.inner.query(&self.file.header, &region)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn data(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data")
            .join(name)
    }

    fn region(contig: &str, start: u64, end: u64) -> ScanRegion {
        ScanRegion {
            contig: contig.to_owned(),
            start,
            end,
        }
    }

    /// Name, flags and 1-based position of the records of a region.
    fn records(file: &AlignmentFile, region: &ScanRegion) -> Vec<(String, u16, usize)> {
        let mut reader = file.reader().unwrap();
        reader
            .query(region)
            .unwrap()
            .map(|record| {
                let record = record.unwrap();
                (
                    record.name().unwrap().to_string(),
                    record.flags().unwrap().bits(),
                    record.alignment_start().unwrap().unwrap().get(),
                )
            })
            .collect()
    }

    #[test]
    fn opens_bam_and_cram() {
        let bam = AlignmentFile::open(&data("positive.bam"), None).unwrap();
        assert_eq!(bam.format(), Format::Bam);
        assert_eq!(bam.index_path(), data("positive.bam.bai"));
        assert_eq!(
            bam.contigs(),
            [
                Contig {
                    name: "chr1".to_owned(),
                    length: 20_000
                },
                Contig {
                    name: "chr2".to_owned(),
                    length: 10_000
                },
            ]
        );
        assert_eq!(bam.samples(), ["positive"]);

        let cram =
            AlignmentFile::open(&data("positive.cram"), Some(&data("reference.fa"))).unwrap();
        assert_eq!(cram.format(), Format::Cram);
        assert_eq!(cram.index_path(), data("positive.cram.crai"));
        assert_eq!(cram.contigs(), bam.contigs());
    }

    #[test]
    fn bam_and_cram_return_the_same_records() {
        let bam = AlignmentFile::open(&data("positive.bam"), None).unwrap();
        let cram =
            AlignmentFile::open(&data("positive.cram"), Some(&data("reference.fa"))).unwrap();
        // Same counts as `samtools view -c` on these regions (unmapped reads
        // placed at their mate included).
        for (region, expected) in [
            (region("chr1", 9_550, 10_450), 703),
            (region("chr1", 2_700, 3_500), 620),
            (region("chr2", 3_700, 4_600), 700),
            (region("chr2", 6_999, 7_281), 153),
        ] {
            let from_bam = records(&bam, &region);
            assert_eq!(from_bam.len(), expected, "{region:?}");
            assert_eq!(from_bam, records(&cram, &region), "{region:?}");
        }
    }

    #[test]
    fn readers_are_independent() {
        let bam = AlignmentFile::open(&data("positive.bam"), None).unwrap();
        let target = region("chr1", 9_550, 10_450);
        let mut first = bam.reader().unwrap();
        let mut second = bam.reader().unwrap();
        let a = first.query(&target).unwrap().count();
        let b = second.query(&target).unwrap().count();
        assert_eq!(a, b);
    }

    #[test]
    fn open_errors() {
        let err = AlignmentFile::open(&data("missing.bam"), None).unwrap_err();
        assert!(matches!(err, AlignmentError::Io { .. }));

        let err = AlignmentFile::open(&data("targets.bed"), None).unwrap_err();
        assert!(matches!(err, AlignmentError::UnknownFormat { .. }));

        let err = AlignmentFile::open(&data("positive.cram"), None).unwrap_err();
        assert_eq!(
            err.to_string(),
            format!(
                "{} is a CRAM file: --reference is required to decode it",
                data("positive.cram").display()
            )
        );

        let err =
            AlignmentFile::open(&data("positive.cram"), Some(&data("missing.fa"))).unwrap_err();
        assert!(matches!(err, AlignmentError::Reference { .. }));
    }

    #[test]
    fn missing_index() {
        let dir = std::env::temp_dir().join(format!("mei-rs-index-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("noindex.bam");
        std::fs::copy(data("positive.bam"), &bam).unwrap();
        let err = AlignmentFile::open(&bam, None).unwrap_err();
        assert!(matches!(err, AlignmentError::MissingIndex { .. }));
        assert!(err.to_string().contains("run `samtools index"));

        // `sample.bai` (without `.bam`) is found too.
        std::fs::copy(data("positive.bam.bai"), dir.join("noindex.bai")).unwrap();
        let file = AlignmentFile::open(&bam, None).unwrap();
        assert_eq!(file.index_path(), dir.join("noindex.bai"));
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
