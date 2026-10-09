//! Local read depth at candidate breakpoints.

use std::io;

use rayon::prelude::*;

use super::alignment::{AlignmentError, AlignmentFile};
use super::filter::ReadFilter;
use crate::regions::targets::ScanRegion;

/// Number of reads passing `filter` that cover one of the two bases around
/// each boundary (`contig`, `boundary`): reads spanning the junction, and
/// reads clipped at it.
///
/// The result is indexed like `boundaries`.
///
/// # Errors
///
/// Fails if the alignment file cannot be read.
pub fn local_depths(
    file: &AlignmentFile,
    boundaries: &[(String, u64)],
    filter: &ReadFilter,
) -> Result<Vec<u64>, AlignmentError> {
    boundaries
        .par_iter()
        .map_init(
            || file.reader(),
            |reader, (contig, boundary)| {
                let reader = reader
                    .as_mut()
                    .map_err(|e| io::Error::other(e.to_string()))?;
                let window = ScanRegion {
                    contig: contig.clone(),
                    start: boundary.saturating_sub(1),
                    end: boundary + 1,
                };
                let mut depth = 0;
                for record in reader.query(&window)? {
                    if filter.check(record?.as_ref())?.is_none() {
                        depth += 1;
                    }
                }
                Ok(depth)
            },
        )
        .collect::<io::Result<Vec<u64>>>()
        .map_err(|source| AlignmentError::Io {
            path: file.path().to_path_buf(),
            source,
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn depth_at_the_junctions() {
        let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let file = AlignmentFile::open(&data.join("positive.bam"), None).unwrap();
        let depths = local_depths(
            &file,
            &[("chr1".to_owned(), 10_000), ("chr2".to_owned(), 4_100)],
            &ReadFilter::default(),
        )
        .unwrap();
        // Same as `samtools view -c -F 0xF04 -q 20` on chr1:10000-10001 and
        // chr2:4100-4101.
        assert_eq!(depths, [155, 121]);
    }
}
