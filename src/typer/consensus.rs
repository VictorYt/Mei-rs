//! MEI consensus sequences and their k-mer index.
//!
//! The default set (`data/consensus/mei_consensus.fa`, built into the binary)
//! holds the human Alu, young LINE-1 and SVA consensus sequences of Dfam
//! (CC0 1.0); `--consensus` replaces it. Headers are `>name#class ...`
//! (`RepeatMasker` style): the class gives the family (`SINE/Alu`, `LINE/L1`,
//! `Retroposon/SVA`); without a class, the name prefix does (`Alu`, `L1`,
//! `SVA`). Sequences of any other family are ignored.

use std::collections::HashMap;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

use thiserror::Error;

/// The default consensus sequences.
pub const DEFAULT_CONSENSUS: &str = include_str!("../../data/consensus/mei_consensus.fa");

/// Length of the k-mers of the index.
pub const K: usize = 11;

/// Mobile element family.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Family {
    /// Alu (SINE).
    Alu,
    /// LINE-1.
    Line1,
    /// SVA.
    Sva,
}

impl Family {
    /// Name, as written in the outputs (and later in the VCF `ALT`,
    /// `<INS:ME:ALU>`).
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::Alu => "ALU",
            Self::Line1 => "LINE1",
            Self::Sva => "SVA",
        }
    }

    /// Family of a FASTA header (`name#class ...`).
    fn of(name: &str, class: Option<&str>) -> Option<Self> {
        let class = class.map(str::to_ascii_lowercase);
        let name = name.to_ascii_lowercase();
        let has = |what: &str| class.as_deref().is_some_and(|c| c.ends_with(what));
        if has("/alu") || (class.is_none() && name.starts_with("alu")) {
            Some(Self::Alu)
        } else if has("/l1") || (class.is_none() && name.starts_with("l1")) {
            Some(Self::Line1)
        } else if has("/sva") || (class.is_none() && name.starts_with("sva")) {
            Some(Self::Sva)
        } else {
            None
        }
    }
}

/// Errors raised while loading consensus sequences.
#[derive(Debug, Error)]
pub enum ConsensusError {
    /// The file cannot be read.
    #[error("cannot read the consensus sequences {path}: {source}")]
    Io {
        /// Path of the FASTA file.
        path: PathBuf,
        /// Underlying I/O error.
        source: io::Error,
    },
    /// No sequence of a known family.
    #[error(
        "no Alu, LINE-1 or SVA sequence in {path} (headers `>name#SINE/Alu`, `>name#LINE/L1`, `>name#Retroposon/SVA`, or names starting with Alu, L1 or SVA)"
    )]
    Empty {
        /// Path of the FASTA file.
        path: PathBuf,
    },
}

/// One consensus sequence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Consensus {
    /// Name (subfamily), e.g. `AluYa5` or `L1HS_3end`.
    pub name: String,
    /// Family.
    pub family: Family,
    /// Bases, upper case.
    pub sequence: Vec<u8>,
}

/// Where a k-mer occurs: consensus index, strand and position (0-based, on
/// the forward consensus for the forward strand, on its reverse complement
/// otherwise).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Occurrence {
    /// Index of the consensus.
    pub consensus: u32,
    /// The k-mer occurs on the reverse complement of the consensus.
    pub reverse: bool,
    /// Start of the k-mer.
    pub position: u32,
}

/// Consensus sequences, with the k-mer index of both strands.
#[derive(Debug)]
pub struct ConsensusSet {
    sequences: Vec<Consensus>,
    reverse: Vec<Vec<u8>>,
    index: HashMap<u32, Vec<Occurrence>>,
}

impl ConsensusSet {
    /// The consensus sequences of `path`, or the default set.
    ///
    /// # Errors
    ///
    /// Fails if the file cannot be read or holds no sequence of a known
    /// family.
    pub fn load(path: Option<&Path>) -> Result<Self, ConsensusError> {
        let Some(path) = path else {
            return Ok(Self::parse(DEFAULT_CONSENSUS));
        };
        let text = fs::read_to_string(path).map_err(|source| ConsensusError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        let set = Self::parse(&text);
        if set.sequences.is_empty() {
            return Err(ConsensusError::Empty {
                path: path.to_path_buf(),
            });
        }
        Ok(set)
    }

    /// Parse FASTA text; sequences of unknown families are skipped.
    #[must_use]
    pub fn parse(text: &str) -> Self {
        let mut sequences = Vec::new();
        let mut current: Option<Consensus> = None;
        for line in text.lines() {
            let line = line.trim_end();
            if let Some(header) = line.strip_prefix('>') {
                sequences.extend(current.take());
                let id = header.split_whitespace().next().unwrap_or("");
                let (name, class) = match id.split_once('#') {
                    Some((name, class)) => (name, Some(class)),
                    None => (id, None),
                };
                current = Family::of(name, class).map(|family| Consensus {
                    name: name.to_owned(),
                    family,
                    sequence: Vec::new(),
                });
            } else if let Some(c) = current.as_mut() {
                c.sequence.extend(
                    line.bytes()
                        .filter(u8::is_ascii_alphabetic)
                        .map(|b| b.to_ascii_uppercase()),
                );
            }
        }
        sequences.extend(current);
        sequences.retain(|c| !c.sequence.is_empty());
        Self::new(sequences)
    }

    fn new(sequences: Vec<Consensus>) -> Self {
        let reverse: Vec<Vec<u8>> = sequences
            .iter()
            .map(|c| reverse_complement(&c.sequence))
            .collect();
        let mut index: HashMap<u32, Vec<Occurrence>> = HashMap::new();
        for (i, c) in sequences.iter().enumerate() {
            for (reverse_strand, sequence) in [(false, &c.sequence), (true, &reverse[i])] {
                for (position, code) in kmers(sequence) {
                    index.entry(code).or_default().push(Occurrence {
                        consensus: u32::try_from(i).unwrap_or(u32::MAX),
                        reverse: reverse_strand,
                        position: u32::try_from(position).unwrap_or(u32::MAX),
                    });
                }
            }
        }
        Self {
            sequences,
            reverse,
            index,
        }
    }

    /// The consensus sequences.
    #[must_use]
    pub fn sequences(&self) -> &[Consensus] {
        &self.sequences
    }

    /// A strand of a consensus: forward, or reverse complement.
    #[must_use]
    pub fn strand(&self, consensus: usize, reverse: bool) -> &[u8] {
        if reverse {
            &self.reverse[consensus]
        } else {
            &self.sequences[consensus].sequence
        }
    }

    /// Occurrences of a k-mer code in the consensus sequences.
    #[must_use]
    pub fn occurrences(&self, code: u32) -> &[Occurrence] {
        self.index.get(&code).map_or(&[], Vec::as_slice)
    }
}

/// Reverse complement (other bases become `N`).
#[must_use]
pub fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

/// The informative k-mers of a sequence: start and 2-bit code. K-mers with
/// a base other than A, C, G or T, or made at least 80% of one base (poly(A)
/// tails, A-rich linkers), are skipped: they seed everywhere.
pub fn kmers(sequence: &[u8]) -> impl Iterator<Item = (usize, u32)> + '_ {
    sequence.windows(K).enumerate().filter_map(|(start, kmer)| {
        let mut code = 0u32;
        let mut counts = [0usize; 4];
        for &b in kmer {
            let bits = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None,
            };
            counts[bits as usize] += 1;
            code = (code << 2) | bits;
        }
        let low_complexity = counts.iter().any(|&n| 5 * n >= 4 * K);
        (!low_complexity).then_some((start, code))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_set() {
        let set = ConsensusSet::load(None).unwrap();
        let count = |family| {
            set.sequences()
                .iter()
                .filter(|c| c.family == family)
                .count()
        };
        assert_eq!(count(Family::Alu), 42);
        assert_eq!(count(Family::Line1), 10);
        assert_eq!(count(Family::Sva), 6);
        let alu = &set.sequences()[0];
        assert_eq!(alu.name, "AluY");
        assert_eq!(alu.sequence.len(), 311);
        assert!(alu.sequence.starts_with(b"GGCCGGGCGCGGTGGCTCAC"));
    }

    #[test]
    fn families_from_headers() {
        let set = ConsensusSet::parse(
            ">AluY#SINE/Alu DF1\nacgt\nACGT\n>MER5#DNA/hAT\nACGT\n>L1HS_3end\nACGT\n>SVA_F#Retroposon/SVA\nAC\nGT\n>x#LINE/L1\n\n",
        );
        let names: Vec<_> = set
            .sequences()
            .iter()
            .map(|c| (c.name.as_str(), c.family, c.sequence.len()))
            .collect();
        assert_eq!(
            names,
            [
                ("AluY", Family::Alu, 8),
                ("L1HS_3end", Family::Line1, 4),
                ("SVA_F", Family::Sva, 4),
            ]
        );
    }

    #[test]
    fn index_both_strands() {
        let set = ConsensusSet::parse(">AluTest#SINE/Alu\nGGCCGGGCGCGGTGGCTCACGCC\n");
        let (_, code) = kmers(b"GGCCGGGCGCG").next().unwrap();
        assert_eq!(
            set.occurrences(code),
            [Occurrence {
                consensus: 0,
                reverse: false,
                position: 0
            }]
        );
        let rc = reverse_complement(b"GGCCGGGCGCG");
        let (_, code) = kmers(&rc).next().unwrap();
        assert_eq!(
            set.occurrences(code),
            [Occurrence {
                consensus: 0,
                reverse: true,
                position: 23 - 11
            }]
        );
    }

    #[test]
    fn low_complexity_kmers_are_skipped() {
        assert_eq!(kmers(b"AAAAAAAAAAA").count(), 0);
        assert_eq!(kmers(b"AAAAAAAAACG").count(), 0); // 9 of 11 A
        assert_eq!(kmers(b"AAAAAAAACGT").count(), 1); // 8 of 11
        assert_eq!(kmers(b"ACGTNACGTACGT").count(), 0);
    }

    #[test]
    fn missing_or_empty_file() {
        let err = ConsensusSet::load(Some(Path::new("missing.fa"))).unwrap_err();
        assert!(matches!(err, ConsensusError::Io { .. }));
        let dir = std::env::temp_dir().join(format!("mei-rs-consensus-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("dna.fa");
        fs::write(&path, ">MER5#DNA/hAT\nACGT\n").unwrap();
        let err = ConsensusSet::load(Some(&path)).unwrap_err();
        assert!(matches!(err, ConsensusError::Empty { .. }));
        fs::remove_dir_all(&dir).unwrap();
    }
}
