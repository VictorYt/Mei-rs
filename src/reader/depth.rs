//! Local read depth at candidate breakpoints, and reads on regions.

use super::alignment::{AlignmentError, AlignmentFile};
use super::batch;
use super::filter::ReadFilter;
use crate::regions::targets::ScanRegion;

/// Number of reads passing `filter` that cover one of the two bases around
/// each boundary (`contig`, `boundary`): reads spanning the junction, and
/// reads clipped at it.
///
/// The result is indexed like `boundaries`. Neighbouring boundaries are read
/// in batches (see [`super::batch`]); each read counts for every boundary it
/// covers.
///
/// # Errors
///
/// Fails if the alignment file cannot be read.
pub fn local_depths(
    file: &AlignmentFile,
    boundaries: &[(String, u64)],
    filter: &ReadFilter,
) -> Result<Vec<u64>, AlignmentError> {
    let windows: Vec<ScanRegion> = boundaries
        .iter()
        .map(|(contig, boundary)| ScanRegion {
            contig: contig.clone(),
            start: boundary.saturating_sub(1),
            end: boundary + 1,
        })
        .collect();
    read_counts(file, &windows, filter)
}

/// Number of reads passing `filter` that overlap each region, indexed like
/// `regions` (in any order). Neighbouring regions are read in batches (see
/// [`super::batch`]); each read counts for every region it overlaps.
///
/// # Errors
///
/// Fails if the alignment file cannot be read.
pub fn read_counts(
    file: &AlignmentFile,
    regions: &[ScanRegion],
    filter: &ReadFilter,
) -> Result<Vec<u64>, AlignmentError> {
    // Sorted by contig and position (contigs in header order are not
    // needed: any grouping works).
    let mut order: Vec<usize> = (0..regions.len()).collect();
    order.sort_by(|&a, &b| {
        let (a, b) = (&regions[a], &regions[b]);
        (&a.contig, a.start, a.end).cmp(&(&b.contig, b.start, b.end))
    });
    let windows: Vec<ScanRegion> = order.iter().map(|&i| regions[i].clone()).collect();
    let counts = batch::fold(
        file,
        &batch::batches(&windows, file.format()),
        |reader, batch| {
            let windows = &windows[batch.range.clone()];
            let mut counts = vec![0; windows.len()];
            for record in reader.query(&batch.span)? {
                let record = record?;
                let covered = batch::overlap(windows, record.as_ref())?;
                if !covered.is_empty() && filter.check(record.as_ref())?.is_none() {
                    for count in &mut counts[covered] {
                        *count += 1;
                    }
                }
            }
            Ok((batch.range.start, counts))
        },
        |mut acc: Vec<(usize, Vec<u64>)>, batch| {
            acc.push(batch);
            acc
        },
        |mut a, b| {
            a.extend(b);
            a
        },
    )?;
    let mut result = vec![0; regions.len()];
    for (first, counts) in counts {
        for (offset, count) in counts.into_iter().enumerate() {
            result[order[first + offset]] = count;
        }
    }
    Ok(result)
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
        assert_eq!(depths, [155, 130]);
    }

    #[test]
    fn cram_batches_give_the_same_depths() {
        let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let bam = AlignmentFile::open(&data.join("positive.bam"), None).unwrap();
        let cram = AlignmentFile::open(
            &data.join("positive.cram"),
            Some(&data.join("reference.fa")),
        )
        .unwrap();
        // Unsorted, close, duplicated and edge boundaries: one CRAM batch per
        // contig, one BAM query each.
        let boundaries: Vec<(String, u64)> = [
            ("chr1", 10_012),
            ("chr2", 4_100),
            ("chr1", 10_000),
            ("chr1", 10_001),
            ("chr1", 3_000),
            ("chr1", 10_000),
            ("chr1", 0),
            ("chr1", 9_400),
        ]
        .into_iter()
        .map(|(c, b)| (c.to_owned(), b))
        .collect();
        let filter = ReadFilter::default();
        let from_bam = local_depths(&bam, &boundaries, &filter).unwrap();
        assert_eq!(from_bam[2], 155);
        assert_eq!(from_bam[2], from_bam[5]);
        assert_eq!(local_depths(&cram, &boundaries, &filter).unwrap(), from_bam);
    }
}
