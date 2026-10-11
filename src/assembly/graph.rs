//! A small de Bruijn graph of the reads of one junction, and the walk that
//! extends a contig from a seed.
//!
//! The graph is the set of k-mers of the reads with their counts (k-mers
//! with a masked base are left out). From a seed (the reference flank next
//! to the junction), the contig is extended one base at a time with the most
//! supported next k-mer; the walk stops when no next k-mer has
//! [`MIN_COUNT`] reads, when two next k-mers are close in support (a branch:
//! another allele or a repeat), on a k-mer already visited (a cycle), or at
//! [`MAX_CONTIG`] bases.
//!
//! A homopolymer longer than k (a long poly(A) tail) is a single k-mer
//! looping on itself, so its length is not in the graph: it is measured on
//! the reads that hold the bases before it, the whole run and the base after
//! it (the median of their run lengths), and the walk resumes after it.

use std::collections::{HashMap, HashSet};

/// Length of the k-mers.
pub const K: usize = 21;

/// Minimum number of reads supporting a k-mer to extend a contig.
pub const MIN_COUNT: u32 = 2;

/// Maximum length of a contig, seed included.
pub const MAX_CONTIG: usize = 600;

/// Bases before a homopolymer that locate it in the reads.
pub const RUN_ANCHOR: usize = 8;

/// 2-bit code of a k-mer, `None` with a base other than `ACGT`.
fn code(kmer: &[u8]) -> Option<u64> {
    kmer.iter().try_fold(0u64, |code, &base| {
        let bits = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        Some(code << 2 | bits)
    })
}

/// k-mers of a set of reads (in the orientation of the walk).
pub struct Graph<'a> {
    counts: HashMap<u64, u32>,
    reads: Vec<&'a [u8]>,
}

impl<'a> Graph<'a> {
    /// Graph of the k-mers of the reads.
    #[must_use]
    pub fn new(reads: Vec<&'a [u8]>) -> Self {
        let mut counts: HashMap<u64, u32> = HashMap::new();
        for read in &reads {
            for kmer in read.windows(K) {
                if let Some(code) = code(kmer) {
                    *counts.entry(code).or_default() += 1;
                }
            }
        }
        Self { counts, reads }
    }

    /// Number of occurrences of a k-mer in the reads.
    #[must_use]
    pub fn count(&self, kmer: &[u8]) -> u32 {
        code(kmer).map_or(0, |c| self.counts.get(&c).copied().unwrap_or(0))
    }

    /// Extend a contig from a seed of [`K`] bases; the contig starts with
    /// the seed.
    #[must_use]
    pub fn extend(&self, seed: &[u8]) -> Vec<u8> {
        let mut contig = seed.to_vec();
        let mut visited: HashSet<u64> = code(seed).into_iter().collect();
        let mut kmer = Vec::with_capacity(K);
        while contig.len() < MAX_CONTIG {
            let prefix = &contig[contig.len() - (K - 1)..];
            let mut next: Vec<(u32, u8)> = b"ACGT"
                .iter()
                .map(|&base| {
                    kmer.clear();
                    kmer.extend_from_slice(prefix);
                    kmer.push(base);
                    (self.count(&kmer), base)
                })
                .collect();
            next.sort_unstable_by(|a, b| b.cmp(a));
            let ((best, base), (second, _)) = (next[0], next[1]);
            if best < MIN_COUNT || (second >= MIN_COUNT && 2 * second >= best) {
                break;
            }
            if prefix.iter().all(|&b| b == base) {
                if !self.extend_run(&mut contig) {
                    break;
                }
                continue;
            }
            kmer.clear();
            kmer.extend_from_slice(prefix);
            kmer.push(base);
            if !code(&kmer).is_some_and(|c| visited.insert(c)) {
                break;
            }
            contig.push(base);
        }
        contig.truncate(MAX_CONTIG);
        contig
    }

    /// Extend a contig ending with a run of `K - 1` identical bases or more
    /// over the whole run and the base after it, measured on the reads;
    /// `false` if fewer than [`MIN_COUNT`] reads hold the whole run, or if
    /// the base after it is masked (`N`, the contig ends there).
    fn extend_run(&self, contig: &mut Vec<u8>) -> bool {
        let base = contig[contig.len() - 1];
        let run = contig.iter().rev().take_while(|&&b| b == base).count();
        let Some(anchor_start) = contig.len().checked_sub(run + RUN_ANCHOR) else {
            return false;
        };
        let anchor = &contig[anchor_start..contig.len() - run];
        // Run length and next base in each read holding the anchor, the run
        // so far and the end of the run.
        let mut found: Vec<(usize, u8)> = Vec::new();
        for read in &self.reads {
            for start in 0..read.len().saturating_sub(RUN_ANCHOR) {
                if &read[start..start + RUN_ANCHOR] != anchor {
                    continue;
                }
                let after = &read[start + RUN_ANCHOR..];
                let length = after.iter().take_while(|&&b| b == base).count();
                if let Some(&next) = after.get(length)
                    && length >= run
                {
                    found.push((length, next));
                }
            }
        }
        if found.len() < MIN_COUNT as usize {
            return false;
        }
        found.sort_unstable();
        let length = found[(found.len() - 1) / 2].0;
        let mut votes: HashMap<u8, usize> = HashMap::new();
        for &(_, next) in &found {
            *votes.entry(next).or_default() += 1;
        }
        let Some((&next, _)) = votes
            .iter()
            .max_by_key(|&(&b, &n)| (n, std::cmp::Reverse(b)))
        else {
            return false;
        };
        contig.extend(std::iter::repeat_n(base, length - run));
        contig.push(next);
        // A masked base ends the run, and the contig.
        next != b'N'
    }
}

/// Reverse complement of a sequence (`N` and other bases kept).
#[must_use]
pub fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reads of `length` bases tiling `sequence` every `step` bases.
    fn tiles(sequence: &[u8], length: usize, step: usize) -> Vec<Vec<u8>> {
        (0..=sequence.len() - length)
            .step_by(step)
            .map(|i| sequence[i..i + length].to_vec())
            .collect()
    }

    const FLANK: &[u8] = b"TGACCTAGGCATTCGAGTCAGCTTAGCCATGGACT";

    /// The contig is the sequence, up to its last base (held by one read
    /// only).
    fn assert_assembled(contig: &[u8], sequence: &[u8]) {
        assert!(
            sequence.starts_with(contig) && contig.len() + 1 >= sequence.len(),
            "{}",
            String::from_utf8_lossy(contig)
        );
    }
    const INSERT: &[u8] = b"GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGG";

    #[test]
    fn extends_into_the_insertion() {
        let alt = [FLANK, INSERT].concat();
        let reads = tiles(&alt, 40, 1);
        let graph = Graph::new(reads.iter().map(Vec::as_slice).collect());
        assert_assembled(&graph.extend(&FLANK[..K]), &alt);
    }

    #[test]
    fn measures_a_long_homopolymer() {
        let tail = [b'A'; 30];
        let alt = [FLANK, &tail[..], b"CGTCAGGATCCATG"].concat();
        let reads = tiles(&alt, 60, 1);
        let graph = Graph::new(reads.iter().map(Vec::as_slice).collect());
        assert_assembled(&graph.extend(&FLANK[..K]), &alt);
    }

    #[test]
    fn stops_without_support_or_at_a_branch() {
        let alt = [FLANK, INSERT].concat();
        // A single read: no k-mer has two reads.
        let single = [alt.clone()];
        let graph = Graph::new(single.iter().map(Vec::as_slice).collect());
        assert_eq!(graph.extend(&FLANK[..K]), &FLANK[..K]);
        // Two alleles after the flank, equally supported.
        let other = [FLANK, b"TTTTCCCCAAAAGGGGTTTTACGTACGT"].concat();
        let mut reads = tiles(&alt, 40, 1);
        reads.extend(tiles(&other, 40, 1));
        let graph = Graph::new(reads.iter().map(Vec::as_slice).collect());
        assert_eq!(graph.extend(&FLANK[..K]), FLANK);
    }

    #[test]
    fn masked_bases_and_reverse_complement() {
        assert_eq!(code(b"ACGN"), None);
        assert_eq!(code(b"ACGT"), Some(0b00_01_10_11));
        assert_eq!(reverse_complement(b"AACGTN"), b"NACGTT");
    }
}
