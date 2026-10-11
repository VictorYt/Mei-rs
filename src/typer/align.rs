//! Local alignment (Smith-Waterman with affine gaps, Gotoh) of a
//! soft-clipped sequence against a window of a consensus sequence.
//!
//! Sequences are short (a clip of at most a read length against a window of
//! a few hundred bases), so a scalar implementation is enough (see the
//! [`crate::typer`] module documentation).

/// Score of a match.
pub const MATCH: i32 = 2;
/// Score of a mismatch (an `N` never matches).
pub const MISMATCH: i32 = -3;
/// Cost of opening a gap: a gap of `k` bases costs `GAP_OPEN + k * GAP_EXTEND`.
pub const GAP_OPEN: i32 = 5;
/// Cost of each gap base.
pub const GAP_EXTEND: i32 = 2;

/// Best local alignment of a query against a target.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Alignment {
    /// Alignment score.
    pub score: i32,
    /// Start of the aligned part of the query, 0-based.
    pub query_start: usize,
    /// End of the aligned part of the query, exclusive.
    pub query_end: usize,
    /// Start of the aligned part of the target, 0-based.
    pub target_start: usize,
    /// End of the aligned part of the target, exclusive.
    pub target_end: usize,
    /// Identical aligned bases.
    pub matches: usize,
    /// Alignment columns: matches, mismatches and gap bases.
    pub columns: usize,
}

impl Alignment {
    /// Fraction of identical columns.
    #[must_use]
    pub fn identity(&self) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        if self.columns == 0 {
            0.0
        } else {
            self.matches as f64 / self.columns as f64
        }
    }
}

const NONE: i32 = i32::MIN / 4;

/// Traceback state: in the H, E or F matrix.
#[derive(Clone, Copy, PartialEq)]
enum State {
    H,
    E,
    F,
}

/// Best local alignment of `query` against `target`, `None` without any
/// positive score. Bases are compared as given (callers pass upper case).
#[must_use]
#[allow(clippy::many_single_char_names)] // the usual names of the recurrence
pub fn local(query: &[u8], target: &[u8]) -> Option<Alignment> {
    let (n, m) = (query.len(), target.len());
    if n == 0 || m == 0 {
        return None;
    }
    let width = m + 1;
    let at = |i: usize, j: usize| i * width + j;
    let score = |i: usize, j: usize| {
        let (q, t) = (query[i - 1], target[j - 1]);
        if q == t && q != b'N' { MATCH } else { MISMATCH }
    };
    // H: best score of an alignment ending at (i, j); E: ending with a gap in
    // the query (consumes the target); F: ending with a gap in the target
    // (consumes the query).
    let mut h = vec![0; (n + 1) * width];
    let mut e = vec![NONE; (n + 1) * width];
    let mut f = vec![NONE; (n + 1) * width];
    let mut best = (0, 0, 0);
    for i in 1..=n {
        for j in 1..=m {
            let cell = at(i, j);
            e[cell] = (h[at(i, j - 1)] - GAP_OPEN - GAP_EXTEND).max(e[at(i, j - 1)] - GAP_EXTEND);
            f[cell] = (h[at(i - 1, j)] - GAP_OPEN - GAP_EXTEND).max(f[at(i - 1, j)] - GAP_EXTEND);
            let value = (h[at(i - 1, j - 1)] + score(i, j))
                .max(e[cell])
                .max(f[cell])
                .max(0);
            h[cell] = value;
            if value > best.0 {
                best = (value, i, j);
            }
        }
    }
    let (best_score, end_i, end_j) = best;
    if best_score <= 0 {
        return None;
    }

    // Traceback, recomputing which move gave each value.
    let (mut i, mut j, mut state) = (end_i, end_j, State::H);
    let (mut matches, mut columns) = (0, 0);
    loop {
        let cell = at(i, j);
        match state {
            State::H => {
                let value = h[cell];
                if value == 0 {
                    break;
                }
                if i > 0 && j > 0 && value == h[at(i - 1, j - 1)] + score(i, j) {
                    if score(i, j) == MATCH {
                        matches += 1;
                    }
                    columns += 1;
                    i -= 1;
                    j -= 1;
                } else if value == e[cell] {
                    state = State::E;
                } else {
                    state = State::F;
                }
            }
            State::E => {
                columns += 1;
                let opened = e[cell] == h[at(i, j - 1)] - GAP_OPEN - GAP_EXTEND;
                j -= 1;
                if opened {
                    state = State::H;
                }
            }
            State::F => {
                columns += 1;
                let opened = f[cell] == h[at(i - 1, j)] - GAP_OPEN - GAP_EXTEND;
                i -= 1;
                if opened {
                    state = State::H;
                }
            }
        }
    }
    Some(Alignment {
        score: best_score,
        query_start: i,
        query_end: end_i,
        target_start: j,
        target_end: end_j,
        matches,
        columns,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_match_inside_the_target() {
        let a = local(b"GATTACA", b"CCCCGATTACACCCC").unwrap();
        assert_eq!(a.score, 7 * MATCH);
        assert_eq!((a.query_start, a.query_end), (0, 7));
        assert_eq!((a.target_start, a.target_end), (4, 11));
        assert_eq!((a.matches, a.columns), (7, 7));
        assert!((a.identity() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn mismatch_and_local_ends() {
        // One mismatch in the middle; the unrelated ends of the query are
        // left out.
        let a = local(b"TTTTGATTGCAGGTACCATTTT", b"AAGATTACAGGTACCAAA").unwrap();
        assert_eq!((a.query_start, a.query_end), (4, 18));
        assert_eq!((a.target_start, a.target_end), (2, 16));
        assert_eq!((a.matches, a.columns), (13, 14));
        assert_eq!(a.score, 13 * MATCH + MISMATCH);
    }

    #[test]
    fn gap() {
        // Three bases deleted from the query.
        let target = b"ACGTTGCAAGGCTTACCGATGCATTGCA";
        let query = b"ACGTTGCAAGGCCCGATGCATTGCA";
        let a = local(query, target).unwrap();
        assert_eq!((a.query_start, a.query_end), (0, query.len()));
        assert_eq!((a.target_start, a.target_end), (0, target.len()));
        assert_eq!(a.matches, query.len());
        assert_eq!(a.columns, target.len());
        assert_eq!(a.score, 25 * MATCH - GAP_OPEN - 3 * GAP_EXTEND);
    }

    #[test]
    fn nothing_to_align() {
        assert_eq!(local(b"", b"ACGT"), None);
        assert_eq!(local(b"AAAA", b"CCCC"), None);
        assert_eq!(local(b"NNNN", b"NNNN"), None);
    }
}
