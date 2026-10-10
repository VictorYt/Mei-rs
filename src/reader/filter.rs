//! Read filtering before signal extraction.
//!
//! Only primary, mapped, unique-enough reads can anchor an insertion signal.
//! Supplementary alignments are not signals of their own: the primary
//! alignment carries them in its `SA` tag.

use std::io;

use noodles::sam::alignment::Record;

/// Read filtering options.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadFilter {
    /// Minimum mapping quality of an anchored read.
    pub min_mapq: u8,
    /// Keep reads flagged as duplicates (for files without duplicate marking).
    pub include_duplicates: bool,
}

impl Default for ReadFilter {
    fn default() -> Self {
        Self {
            min_mapq: 20,
            include_duplicates: false,
        }
    }
}

/// Why a read was set aside.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Exclusion {
    /// Unmapped (it may still be the mate of an anchored read).
    Unmapped,
    /// Secondary alignment.
    Secondary,
    /// Supplementary alignment.
    Supplementary,
    /// Failed the vendor quality checks.
    QcFail,
    /// PCR or optical duplicate.
    Duplicate,
    /// Mapping quality below the threshold.
    LowMapq,
}

impl ReadFilter {
    /// `None` if the read can anchor a signal, or the reason it cannot.
    ///
    /// # Errors
    ///
    /// Fails if the flags or the mapping quality cannot be decoded.
    pub fn check(&self, record: &dyn Record) -> io::Result<Option<Exclusion>> {
        let flags = record.flags()?;
        let exclusion = if flags.is_unmapped() {
            Some(Exclusion::Unmapped)
        } else if flags.is_secondary() {
            Some(Exclusion::Secondary)
        } else if flags.is_supplementary() {
            Some(Exclusion::Supplementary)
        } else if flags.is_qc_fail() {
            Some(Exclusion::QcFail)
        } else if flags.is_duplicate() && !self.include_duplicates {
            Some(Exclusion::Duplicate)
        } else if mapq(record)?.is_some_and(|q| q < self.min_mapq) {
            Some(Exclusion::LowMapq)
        } else {
            None
        };
        Ok(exclusion)
    }
}

/// Mapping quality, `None` if unavailable (255).
///
/// # Errors
///
/// Fails if the mapping quality cannot be decoded.
pub fn mapq(record: &dyn Record) -> io::Result<Option<u8>> {
    record
        .mapping_quality()
        .transpose()
        .map(|q| q.map(u8::from))
}

/// Read counts, by outcome of the filter.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FilterStats {
    /// Reads read.
    pub total: u64,
    /// Reads kept as possible anchors.
    pub kept: u64,
    /// Unmapped reads.
    pub unmapped: u64,
    /// Secondary alignments.
    pub secondary: u64,
    /// Supplementary alignments.
    pub supplementary: u64,
    /// Reads failing the vendor quality checks.
    pub qc_fail: u64,
    /// Duplicates.
    pub duplicate: u64,
    /// Reads below the mapping quality threshold.
    pub low_mapq: u64,
}

impl FilterStats {
    /// Count a read.
    pub fn add(&mut self, exclusion: Option<Exclusion>) {
        self.total += 1;
        let counter = match exclusion {
            None => &mut self.kept,
            Some(Exclusion::Unmapped) => &mut self.unmapped,
            Some(Exclusion::Secondary) => &mut self.secondary,
            Some(Exclusion::Supplementary) => &mut self.supplementary,
            Some(Exclusion::QcFail) => &mut self.qc_fail,
            Some(Exclusion::Duplicate) => &mut self.duplicate,
            Some(Exclusion::LowMapq) => &mut self.low_mapq,
        };
        *counter += 1;
    }

    /// Add the counts of `other`.
    pub fn merge(&mut self, other: &Self) {
        self.total += other.total;
        self.kept += other.kept;
        self.unmapped += other.unmapped;
        self.secondary += other.secondary;
        self.supplementary += other.supplementary;
        self.qc_fail += other.qc_fail;
        self.duplicate += other.duplicate;
        self.low_mapq += other.low_mapq;
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use noodles::core::Position;
    use noodles::sam::alignment::RecordBuf;
    use noodles::sam::alignment::record::{Flags, MappingQuality};

    pub(crate) fn record(flags: Flags, mapq: u8) -> RecordBuf {
        let mut builder = RecordBuf::builder()
            .set_name("r1")
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::MIN);
        if let Some(mapq) = MappingQuality::new(mapq) {
            builder = builder.set_mapping_quality(mapq);
        }
        builder.build()
    }

    #[test]
    fn exclusions() {
        let filter = ReadFilter::default();
        let check = |flags, mapq| filter.check(&record(flags, mapq)).unwrap();
        assert_eq!(check(Flags::SEGMENTED, 60), None);
        assert_eq!(check(Flags::UNMAPPED, 60), Some(Exclusion::Unmapped));
        assert_eq!(check(Flags::SECONDARY, 60), Some(Exclusion::Secondary));
        assert_eq!(
            check(Flags::SUPPLEMENTARY, 60),
            Some(Exclusion::Supplementary)
        );
        assert_eq!(check(Flags::QC_FAIL, 60), Some(Exclusion::QcFail));
        assert_eq!(check(Flags::DUPLICATE, 60), Some(Exclusion::Duplicate));
        assert_eq!(check(Flags::empty(), 19), Some(Exclusion::LowMapq));
        assert_eq!(check(Flags::empty(), 20), None);
        // MAPQ 255 means "unavailable": kept.
        assert_eq!(check(Flags::empty(), 255), None);
        // Unmapped takes precedence.
        assert_eq!(
            check(Flags::UNMAPPED | Flags::DUPLICATE, 0),
            Some(Exclusion::Unmapped)
        );
    }

    #[test]
    fn include_duplicates() {
        let filter = ReadFilter {
            include_duplicates: true,
            ..ReadFilter::default()
        };
        assert_eq!(filter.check(&record(Flags::DUPLICATE, 60)).unwrap(), None);
        assert_eq!(
            filter.check(&record(Flags::DUPLICATE, 0)).unwrap(),
            Some(Exclusion::LowMapq)
        );
    }

    #[test]
    fn stats() {
        let mut stats = FilterStats::default();
        stats.add(None);
        stats.add(Some(Exclusion::Duplicate));
        stats.add(Some(Exclusion::Duplicate));
        let mut total = FilterStats::default();
        total.merge(&stats);
        total.merge(&stats);
        assert_eq!(total.total, 6);
        assert_eq!(total.kept, 2);
        assert_eq!(total.duplicate, 4);
    }
}
