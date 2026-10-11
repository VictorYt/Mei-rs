//! Unmapped mates of the candidates, read again after typing.
//!
//! An unmapped mate is placed at its anchor (same contig and position), so
//! the records at the anchors of a candidate hold its unmapped mates. They
//! are read once the candidates are typed, for the mates to check only (see
//! [`crate::typer::mates`]), instead of keeping every unmapped read during
//! the extraction.
//!
//! Coordinates: [`MateRequest::anchor_start`] is 0-based; a placed unmapped
//! record has the 1-based position of its anchor.

use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasher, BuildHasherDefault, DefaultHasher};

use super::alignment::{AlignmentError, AlignmentFile};
use super::batch;
use crate::regions::targets::ScanRegion;

/// Hash of a read name, to pair unmapped mates with their anchors.
#[must_use]
pub fn name_hash(name: &[u8]) -> u64 {
    BuildHasherDefault::<DefaultHasher>::default().hash_one(name)
}

/// An unmapped mate to read: its anchor.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MateRequest {
    /// Contig of the anchor (index in the header).
    pub contig: usize,
    /// 0-based alignment start of the anchor (and position of the mate).
    pub anchor_start: u64,
    /// Hash of the read name ([`name_hash`]).
    pub name: u64,
}

/// Sequences of the unmapped mates requested, by group (one group per
/// candidate, all on one contig), in request order: `None` when not found.
///
/// # Errors
///
/// Fails if the alignment file cannot be read.
pub fn fetch_unmapped_mates(
    file: &AlignmentFile,
    contigs: &[String],
    groups: &[Vec<MateRequest>],
) -> Result<Vec<Vec<Option<Vec<u8>>>>, AlignmentError> {
    // One window per group: from its first to its last anchor start.
    let mut order: Vec<(usize, ScanRegion)> = groups
        .iter()
        .enumerate()
        .filter_map(|(g, requests)| {
            let first = requests.iter().map(|r| r.anchor_start).min()?;
            let last = requests.iter().map(|r| r.anchor_start).max()?;
            Some((
                g,
                ScanRegion {
                    contig: contigs[requests[0].contig].clone(),
                    start: first,
                    end: last + 1,
                },
            ))
        })
        .collect();
    order.sort_by(|a, b| (&a.1.contig, a.1.start).cmp(&(&b.1.contig, b.1.start)));
    let windows: Vec<ScanRegion> = order.iter().map(|(_, w)| w.clone()).collect();
    let names: Vec<HashSet<u64>> = order
        .iter()
        .map(|&(g, _)| groups[g].iter().map(|r| r.name).collect())
        .collect();

    let found = batch::fold(
        file,
        &batch::batches(&windows, file.format()),
        |reader, batch| {
            let mut found: Vec<(usize, u64, Vec<u8>)> = Vec::new();
            for record in reader.query(&batch.span)? {
                let record = record?;
                let flags = record.flags()?;
                if !flags.is_unmapped() || flags.is_mate_unmapped() {
                    continue;
                }
                let (Some(name), Some(start)) =
                    (record.name(), record.alignment_start().transpose()?)
                else {
                    continue;
                };
                let start = start.get() as u64 - 1;
                let hash = name_hash(name.as_ref());
                for w in batch.range.clone() {
                    let window = &windows[w];
                    if (window.start..window.end).contains(&start) && names[w].contains(&hash) {
                        let mut sequence: Vec<u8> = record.sequence().iter().collect();
                        sequence.make_ascii_uppercase();
                        found.push((w, hash, sequence));
                    }
                }
            }
            Ok(found)
        },
        |mut acc: Vec<(usize, u64, Vec<u8>)>, batch| {
            acc.extend(batch);
            acc
        },
        |mut a, b| {
            a.extend(b);
            a
        },
    )?;

    let mut by_window: Vec<HashMap<u64, Vec<u8>>> = vec![HashMap::new(); windows.len()];
    for (w, hash, sequence) in found {
        by_window[w].insert(hash, sequence);
    }
    let mut sequences: Vec<Vec<Option<Vec<u8>>>> =
        groups.iter().map(|g| vec![None; g.len()]).collect();
    for (w, &(g, _)) in order.iter().enumerate() {
        for (request, sequence) in groups[g].iter().zip(&mut sequences[g]) {
            *sequence = by_window[w].remove(&request.name);
        }
    }
    Ok(sequences)
}
