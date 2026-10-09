//! Parsing of the capture BED file.
//!
//! At least 3 columns are required (contig, start, end); the 4th, when
//! present, is the probe name and the others are ignored. Columns are
//! separated by tabs or spaces. Empty lines, comments (`#`) and `track` /
//! `browser` lines are skipped. Gzip and bgzip compressed files are accepted.

use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

use flate2::bufread::MultiGzDecoder;
use thiserror::Error;

/// A capture probe (one BED line), 0-based half-open.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Probe {
    /// Contig name.
    pub contig: String,
    /// First base, 0-based.
    pub start: u64,
    /// Base after the last one, 0-based.
    pub end: u64,
    /// Name (4th column), if any.
    pub name: Option<String>,
}

/// Errors raised while reading a BED file.
#[derive(Debug, Error)]
pub enum BedError {
    /// The file could not be read.
    #[error("cannot read BED file {path}: {source}")]
    Io {
        /// Path of the BED file.
        path: PathBuf,
        /// Underlying I/O error.
        source: io::Error,
    },
    /// A line has fewer than 3 columns.
    #[error("{path}, line {line}: expected at least 3 columns (contig, start, end), found {found}")]
    MissingColumns {
        /// Path of the BED file.
        path: PathBuf,
        /// 1-based line number.
        line: usize,
        /// Number of columns found.
        found: usize,
    },
    /// A coordinate is not a non-negative integer.
    #[error("{path}, line {line}: invalid {field} `{value}`, expected a non-negative integer")]
    InvalidCoordinate {
        /// Path of the BED file.
        path: PathBuf,
        /// 1-based line number.
        line: usize,
        /// Column name (`start` or `end`).
        field: &'static str,
        /// Value found.
        value: String,
    },
    /// The interval is empty or reversed.
    #[error("{path}, line {line}: start ({start}) must be lower than end ({end})")]
    EmptyInterval {
        /// Path of the BED file.
        path: PathBuf,
        /// 1-based line number.
        line: usize,
        /// Start found.
        start: u64,
        /// End found.
        end: u64,
    },
    /// The file holds no probe.
    #[error("BED file {path} holds no target")]
    Empty {
        /// Path of the BED file.
        path: PathBuf,
    },
}

/// Read the probes of a BED file, in file order.
///
/// # Errors
///
/// Fails if the file cannot be read, if a line is malformed (the error gives
/// its number), or if the file holds no probe.
pub fn read(path: &Path) -> Result<Vec<Probe>, BedError> {
    let io_error = |source| BedError::Io {
        path: path.to_path_buf(),
        source,
    };
    let mut reader = File::open(path).map(BufReader::new).map_err(io_error)?;
    let gzipped = reader
        .fill_buf()
        .map_err(io_error)?
        .starts_with(&[0x1f, 0x8b]);
    if gzipped {
        parse(BufReader::new(MultiGzDecoder::new(reader)), path)
    } else {
        parse(reader, path)
    }
}

/// Parse BED content; `path` is only used in error messages.
///
/// # Errors
///
/// See [`read`].
pub fn parse<R: BufRead>(reader: R, path: &Path) -> Result<Vec<Probe>, BedError> {
    let mut probes = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line.map_err(|source| BedError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        if let Some(probe) = parse_line(&line, i + 1, path)? {
            probes.push(probe);
        }
    }
    if probes.is_empty() {
        return Err(BedError::Empty {
            path: path.to_path_buf(),
        });
    }
    Ok(probes)
}

fn parse_line(line: &str, number: usize, path: &Path) -> Result<Option<Probe>, BedError> {
    let line = line.trim();
    if line.is_empty()
        || line.starts_with('#')
        || line.starts_with("track")
        || line.starts_with("browser")
    {
        return Ok(None);
    }
    let fields: Vec<&str> = line.split_ascii_whitespace().collect();
    if fields.len() < 3 {
        return Err(BedError::MissingColumns {
            path: path.to_path_buf(),
            line: number,
            found: fields.len(),
        });
    }
    let coordinate = |field: &'static str, value: &str| {
        value
            .parse::<u64>()
            .map_err(|_| BedError::InvalidCoordinate {
                path: path.to_path_buf(),
                line: number,
                field,
                value: value.to_owned(),
            })
    };
    let start = coordinate("start", fields[1])?;
    let end = coordinate("end", fields[2])?;
    if start >= end {
        return Err(BedError::EmptyInterval {
            path: path.to_path_buf(),
            line: number,
            start,
            end,
        });
    }
    Ok(Some(Probe {
        contig: fields[0].to_owned(),
        start,
        end,
        name: fields.get(3).map(|&name| name.to_owned()),
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn parse_str(content: &str) -> Result<Vec<Probe>, BedError> {
        parse(content.as_bytes(), Path::new("t.bed"))
    }

    fn probe(contig: &str, start: u64, end: u64, name: Option<&str>) -> Probe {
        Probe {
            contig: contig.to_owned(),
            start,
            end,
            name: name.map(str::to_owned),
        }
    }

    #[test]
    fn parses_probes_and_skips_headers() {
        let content = "track name=panel\nbrowser position chr1\n# comment\n\n\
                       chr1\t100\t200\tEXON1\t0\t+\nchr2 5 10\r\n";
        assert_eq!(
            parse_str(content).unwrap(),
            [
                probe("chr1", 100, 200, Some("EXON1")),
                probe("chr2", 5, 10, None)
            ]
        );
    }

    #[test]
    fn errors_give_the_line_number() {
        let err = parse_str("chr1\t1\t2\nchr1\t5\n").unwrap_err();
        assert_eq!(
            err.to_string(),
            "t.bed, line 2: expected at least 3 columns (contig, start, end), found 2"
        );
        let err = parse_str("chr1\t-1\t2\n").unwrap_err();
        assert_eq!(
            err.to_string(),
            "t.bed, line 1: invalid start `-1`, expected a non-negative integer"
        );
        let err = parse_str("chr1\t1\tx\n").unwrap_err();
        assert!(matches!(
            err,
            BedError::InvalidCoordinate { field: "end", .. }
        ));
        let err = parse_str("chr1\t10\t10\n").unwrap_err();
        assert_eq!(
            err.to_string(),
            "t.bed, line 1: start (10) must be lower than end (10)"
        );
    }

    #[test]
    fn empty_bed_is_an_error() {
        let err = parse_str("# only a comment\n").unwrap_err();
        assert_eq!(err.to_string(), "BED file t.bed holds no target");
    }

    #[test]
    fn reads_plain_and_gzipped_files() {
        let dir = std::env::temp_dir().join(format!("mei-rs-bed-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let content = "chr1\t100\t200\tA\n";

        let plain = dir.join("t.bed");
        std::fs::write(&plain, content).unwrap();
        let gz = dir.join("t.bed.gz");
        let mut encoder =
            flate2::write::GzEncoder::new(File::create(&gz).unwrap(), flate2::Compression::fast());
        encoder.write_all(content.as_bytes()).unwrap();
        encoder.finish().unwrap();

        let expected = [probe("chr1", 100, 200, Some("A"))];
        assert_eq!(read(&plain).unwrap(), expected);
        assert_eq!(read(&gz).unwrap(), expected);
        assert!(matches!(
            read(&dir.join("missing.bed")),
            Err(BedError::Io { .. })
        ));
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
