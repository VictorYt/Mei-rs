//! Genomic regions: the `--region` option, and later the capture targets
//! (BED loading, padding, interval tree).
//!
//! Coordinates are stored 0-based half-open, like BED; they are read and
//! displayed 1-based inclusive, like samtools regions.

use std::fmt;
use std::str::FromStr;

/// A contig, or an interval of a contig.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Region {
    /// Contig name, as in the alignment header.
    pub contig: String,
    /// Interval of the contig; `None` for the whole contig.
    pub span: Option<Span>,
}

/// A 0-based half-open interval (`start < end`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Span {
    /// First base, 0-based.
    pub start: u64,
    /// Base after the last one, 0-based.
    pub end: u64,
}

/// Error returned when a region string cannot be parsed.
#[derive(Debug, PartialEq, Eq)]
pub struct ParseRegionError(String);

impl fmt::Display for ParseRegionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ParseRegionError {}

impl FromStr for Region {
    type Err = ParseRegionError;

    /// Parse `chr` or `chr:start-end` (1-based inclusive, commas allowed in
    /// numbers). A suffix without `-` is part of the contig name, so names
    /// containing `:` (e.g. `HLA-A*01:01`) still work.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (contig, span) = match s.rsplit_once(':') {
            Some((contig, range)) if looks_like_range(range) => (contig, Some(parse_span(range)?)),
            _ => (s, None),
        };
        if contig.is_empty() {
            return Err(ParseRegionError(format!("missing contig name in {s:?}")));
        }
        Ok(Self {
            contig: contig.to_owned(),
            span,
        })
    }
}

impl fmt::Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.span {
            Some(span) => write!(f, "{}:{}-{}", self.contig, span.start + 1, span.end),
            None => f.write_str(&self.contig),
        }
    }
}

fn looks_like_range(range: &str) -> bool {
    range.contains('-')
        && range
            .chars()
            .all(|c| c.is_ascii_digit() || c == ',' || c == '-')
}

fn parse_span(range: &str) -> Result<Span, ParseRegionError> {
    let invalid = || ParseRegionError(format!("invalid range {range:?}, expected start-end"));
    let (start, end) = range.split_once('-').ok_or_else(invalid)?;
    let number = |s: &str| s.replace(',', "").parse::<u64>().map_err(|_| invalid());
    let (start, end) = (number(start)?, number(end)?);
    if start == 0 {
        return Err(ParseRegionError(format!(
            "invalid range {range:?}: positions are 1-based"
        )));
    }
    if end < start {
        return Err(ParseRegionError(format!(
            "invalid range {range:?}: end is before start"
        )));
    }
    Ok(Span {
        start: start - 1,
        end,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region(contig: &str, span: Option<(u64, u64)>) -> Region {
        Region {
            contig: contig.to_owned(),
            span: span.map(|(start, end)| Span { start, end }),
        }
    }

    #[test]
    fn parses_whole_contig() {
        assert_eq!("chr1".parse(), Ok(region("chr1", None)));
    }

    #[test]
    fn parses_interval_as_zero_based_half_open() {
        assert_eq!("chr1:10-20".parse(), Ok(region("chr1", Some((9, 20)))));
        assert_eq!("chr1:5-5".parse(), Ok(region("chr1", Some((4, 5)))));
        assert_eq!(
            "chrX:1,000-2,000".parse(),
            Ok(region("chrX", Some((999, 2000))))
        );
    }

    #[test]
    fn colon_in_contig_name() {
        assert_eq!("HLA-A*01:01".parse(), Ok(region("HLA-A*01:01", None)));
        // Without `-`, the suffix is part of the name; an unknown contig is
        // reported once the alignment header is read.
        assert_eq!("chr1:10".parse(), Ok(region("chr1:10", None)));
        assert_eq!(
            "HLA-A*01:01:1-100".parse(),
            Ok(region("HLA-A*01:01", Some((0, 100))))
        );
    }

    #[test]
    fn rejects_invalid_regions() {
        for s in [
            "",
            ":1-10",
            "chr1:0-10",
            "chr1:20-10",
            "chr1:1-2-3",
            "chr1:-5",
        ] {
            assert!(s.parse::<Region>().is_err(), "{s:?} should be rejected");
        }
    }

    #[test]
    fn display_round_trips() {
        for s in ["chr1", "chr1:10-20", "HLA-A*01:01:1-100"] {
            assert_eq!(s.parse::<Region>().unwrap().to_string(), s);
        }
    }
}
