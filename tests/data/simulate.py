"""Simulate the synthetic test fixtures (see README.md).

Writes, into the directory given as the only argument:
- reference.fa: two random contigs; chr2 holds a reference copy of the Alu;
- targets.bed: three capture targets;
- targets_tiled.bed: the same, with the insertion target split into two
  overlapping probes whose edges fall on the insertion junctions;
- positive.sam: heterozygous Alu insertion with a TSD and a poly(A) tail on
  chr1, a probe-edge artefact and a few marked PCR duplicates;
- negative.sam: same targets, without the insertion;
- truth.tsv: the simulated insertion.

Reads are paired-end 2x100 bp and are "aligned" from the simulation truth:
anchored parts are matches, parts inside the insertion are soft-clipped, reads
entirely inside the Alu map to the chr2 copy with MAPQ 0, and the other reads
inside the insertion are unmapped. Only `random.random()` is used, so the
output does not depend on the Python version.
"""

import math
import os
import random
import sys

SEED = 20261009
READ_LEN = 100
MIN_ANCHOR = 20
DEPTH_PER_HAPLOTYPE = 50
ERROR_RATE = 0.002

# Approximate AluY consensus, for test purposes only: real consensus sequences
# (Dfam) are used for typing.
ALU = (
    "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTC"
    "AGGAGATCGAGACCATCCCGGCTAAAACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGG"
    "GCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGA"
    "GGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTC"
)
POLYA = "A" * 30
INSERTION = ALU + POLYA

CONTIGS = [("chr1", 20000), ("chr2", 10000)]
# Insertion on chr1, after the 0-based position INS_POS (plus strand). The
# target site chr1[INS_POS:INS_POS + TSD_LEN] is duplicated:
# alt = chr1[:INS_POS + TSD_LEN] + INSERTION + chr1[INS_POS:].
INS_POS = 10000
TSD_LEN = 12
# Reference Alu copy on chr2, outside every target.
ALU_COPY_POS = 7000
# Probe whose start produces artefactual soft-clips.
PROBE_EDGE = ("chr1", 3000)
TARGETS = [
    ("chr1", 9850, 10150, "target_insertion"),
    ("chr1", 3000, 3200, "target_probe_edge"),
    ("chr2", 4000, 4300, "target_negative"),
]
# Same targets, the insertion one tiled with two overlapping probes: the first
# ends at the left junction, the second starts at the right one. Every
# soft-clip of the insertion stops at a probe edge, yet its poly(A) tail and
# TSD must keep it from the probe-edge filter.
TILED_TARGETS = [
    ("chr1", 9850, INS_POS + TSD_LEN, "target_insertion_left"),
    ("chr1", INS_POS, 10150, "target_insertion_right"),
] + TARGETS[1:]

COMPLEMENT = str.maketrans("ACGT", "TGCA")


class Rng:
    """Deterministic helpers built on `random.random()` only."""

    def __init__(self, seed):
        self.rng = random.Random(seed)

    def random(self):
        return self.rng.random()

    def below(self, n):
        return min(int(self.rng.random() * n), n - 1)

    def base(self):
        return "ACGT"[self.below(4)]

    def sequence(self, n):
        return "".join(self.base() for _ in range(n))

    def gauss(self, mu, sigma):
        u1 = 1.0 - self.rng.random()
        u2 = self.rng.random()
        return mu + sigma * math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


class Mapping:
    """Alignment of one read: `pos` is 0-based, `None` when unmapped."""

    def __init__(self, contig=None, pos=None, cigar="*", mapq=0, reverse=False):
        self.contig = contig
        self.pos = pos
        self.cigar = cigar
        self.mapq = mapq
        self.reverse = reverse

    @property
    def mapped(self):
        return self.contig is not None

    @property
    def end(self):
        """0-based end on the reference (exclusive)."""
        matched = 0
        number = ""
        for c in self.cigar:
            if c.isdigit():
                number += c
            else:
                if c == "M":
                    matched += int(number)
                number = ""
        return self.pos + matched


class Haplotype:
    """A haplotype sequence of one contig, and how to map it back."""

    def __init__(self, contig, sequence, insertion_at=None):
        self.contig = contig
        self.sequence = sequence
        # alt index where the insertion starts (None for a reference haplotype)
        self.insertion_at = insertion_at

    def to_alt(self, ref_pos):
        """Haplotype index of a reference position (second TSD copy side)."""
        if self.insertion_at is None or ref_pos < INS_POS + TSD_LEN:
            return ref_pos
        return ref_pos + TSD_LEN + len(INSERTION)

    def align(self, start, end):
        """Map the haplotype interval [start, end) to the reference."""
        if self.insertion_at is None:
            return Mapping(self.contig, start, f"{end - start}M", 60)
        ins_start = self.insertion_at
        ins_end = ins_start + len(INSERTION)
        shift = len(INSERTION) + TSD_LEN
        if end <= ins_start:
            return Mapping(self.contig, start, f"{end - start}M", 60)
        if start >= ins_end:
            return Mapping(self.contig, start - shift, f"{end - start}M", 60)
        left = ins_start - start
        right = end - ins_end
        if left >= MIN_ANCHOR:
            return Mapping(self.contig, start, f"{left}M{end - start - left}S", 60)
        if right >= MIN_ANCHOR:
            return Mapping(self.contig, INS_POS, f"{end - start - right}S{right}M", 60)
        if start >= ins_start and end <= ins_start + len(ALU):
            # Multi-mapping read: placed on the reference copy, MAPQ 0.
            offset = start - ins_start
            return Mapping("chr2", ALU_COPY_POS + offset, f"{end - start}M", 0)
        return Mapping()


def sequencing_errors(rng, seq):
    out = []
    for b in seq:
        if rng.random() < ERROR_RATE:
            b = "ACGT".replace(b, "")[rng.below(3)]
        out.append(b)
    return "".join(out)


def qualities(rng, n):
    return "".join(chr(33 + 30 + rng.below(11)) for _ in range(n))


class Writer:
    def __init__(self, path, sample, rng):
        self.out = open(path, "w")
        self.sample = sample
        self.rng = rng
        self.count = 0
        self.out.write("@HD\tVN:1.6\tSO:unsorted\n")
        for name, length in CONTIGS:
            self.out.write(f"@SQ\tSN:{name}\tLN:{length}\n")
        self.out.write(f"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA\n")

    def close(self):
        self.out.close()

    def pair(self, left, right, duplicate=False):
        """Write a fragment. `left` and `right` are (mapping, sequence) on the
        forward strand of the haplotype; `right` is sequenced in reverse."""
        (lmap, lseq), (rmap, rseq) = left, right
        if not lmap.mapped and not rmap.mapped:
            return
        rmap.reverse = True
        self.count += 1
        name = f"{self.sample}_{self.count:05d}"
        forward_first = self.rng.random() < 0.5
        reads = [(lmap, lseq), (rmap, rseq)]
        quals = [qualities(self.rng, len(lseq)), qualities(self.rng, len(rseq))]
        for i, ((m, seq), qual) in enumerate(zip(reads, quals)):
            mate = reads[1 - i][0]
            first = (i == 0) == forward_first
            self.record(name, m, mate, seq, qual, first, duplicate)

    def record(self, name, m, mate, seq, qual, first, duplicate):
        flag = 0x1 | (0x40 if first else 0x80)
        if duplicate:
            flag |= 0x400
        if m.mapped:
            if m.reverse:
                flag |= 0x10
        else:
            flag |= 0x4
            if m.reverse:
                # Unmapped reads keep the sequenced orientation.
                seq, qual = revcomp(seq), qual[::-1]
        if mate.mapped:
            if mate.reverse:
                flag |= 0x20
        else:
            flag |= 0x8
        both = m.mapped and mate.mapped
        proper = (
            both
            and m.contig == mate.contig
            and m.mapq > 0
            and mate.mapq > 0
            and abs(self.tlen(m, mate)) < 1000
        )
        if proper:
            flag |= 0x2
        # Unmapped reads are placed at their mate.
        here = m if m.mapped else mate
        there = mate if mate.mapped else m
        rname, pos = here.contig, here.pos + 1
        rnext = "=" if there.contig == rname else there.contig
        pnext = there.pos + 1
        tlen = self.tlen(m, mate) if both and m.contig == mate.contig else 0
        mapq = m.mapq if m.mapped else 0
        cigar = m.cigar if m.mapped else "*"
        tags = [f"RG:Z:{self.sample}"]
        if mate.mapped:
            tags.append(f"MQ:i:{mate.mapq}")
        fields = [name, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual]
        self.out.write("\t".join(map(str, fields + tags)) + "\n")

    @staticmethod
    def tlen(m, mate):
        start = min(m.pos, mate.pos)
        end = max(m.end, mate.end)
        size = end - start
        leftmost = m.pos < mate.pos or (m.pos == mate.pos and not m.reverse)
        return size if leftmost else -size


def fragments(rng, haplotype, target, writer, duplicates=0):
    """Capture fragments overlapping the target on one haplotype."""
    contig, tstart, tend, _ = target
    start = max(0, haplotype.to_alt(tstart) - 350)
    end = min(len(haplotype.sequence), haplotype.to_alt(tend) + 50)
    count = DEPTH_PER_HAPLOTYPE * (end - start) // (2 * READ_LEN)
    junction_fragments = []
    for _ in range(count):
        length = max(200, min(450, int(round(rng.gauss(300, 30)))))
        fstart = start + rng.below(max(1, end - start - length))
        fend = min(fstart + length, len(haplotype.sequence))
        lstart, lend = fstart, fstart + READ_LEN
        rstart, rend = fend - READ_LEN, fend
        left = (haplotype.align(lstart, lend), sequencing_errors(rng, haplotype.sequence[lstart:lend]))
        right = (haplotype.align(rstart, rend), sequencing_errors(rng, haplotype.sequence[rstart:rend]))
        writer.pair(left, right)
        if "S" in left[0].cigar or "S" in right[0].cigar:
            junction_fragments.append((left, right))
    for left, right in junction_fragments[:duplicates]:
        copy = lambda read: (Mapping(read[0].contig, read[0].pos, read[0].cigar, read[0].mapq), read[1])
        writer.pair(copy(left), copy(right), duplicate=True)


def probe_edge_artefact(rng, reference, writer, count=10):
    """Reads starting exactly at the probe start, with a chimeric soft-clip."""
    contig, pos = PROBE_EDGE
    seq = reference[contig]
    clip = 25
    for _ in range(count):
        length = max(200, min(450, int(round(rng.gauss(300, 30)))))
        anchored = seq[pos : pos + READ_LEN - clip]
        left = (Mapping(contig, pos, f"{clip}S{READ_LEN - clip}M", 60), rng.sequence(clip) + anchored)
        rstart = pos + length - READ_LEN
        right = (Mapping(contig, rstart, f"{READ_LEN}M", 60), sequencing_errors(rng, seq[rstart : rstart + READ_LEN]))
        writer.pair(left, right)


def main(out_dir):
    rng = Rng(SEED)
    reference = {name: rng.sequence(length) for name, length in CONTIGS}
    chr2 = reference["chr2"]
    reference["chr2"] = chr2[:ALU_COPY_POS] + ALU + chr2[ALU_COPY_POS + len(ALU) :]

    with open(os.path.join(out_dir, "reference.fa"), "w") as fa:
        for name, _ in CONTIGS:
            fa.write(f">{name}\n")
            seq = reference[name]
            for i in range(0, len(seq), 60):
                fa.write(seq[i : i + 60] + "\n")

    for file_name, targets in (("targets.bed", TARGETS), ("targets_tiled.bed", TILED_TARGETS)):
        with open(os.path.join(out_dir, file_name), "w") as bed:
            for contig, start, end, name in targets:
                bed.write(f"{contig}\t{start}\t{end}\t{name}\n")

    chr1 = reference["chr1"]
    alt = chr1[: INS_POS + TSD_LEN] + INSERTION + chr1[INS_POS:]
    tsd = chr1[INS_POS : INS_POS + TSD_LEN]
    with open(os.path.join(out_dir, "truth.tsv"), "w") as truth:
        truth.write("#contig\tleft_junction\tright_junction\ttsd_start\ttsd_end\ttsd\telement\tlength\tpolya\tstrand\tgenotype\n")
        truth.write(
            f"chr1\t{INS_POS + TSD_LEN}\t{INS_POS + 1}\t{INS_POS + 1}\t{INS_POS + TSD_LEN}\t{tsd}"
            f"\tAluY\t{len(ALU)}\t{len(POLYA)}\t+\t0/1\n"
        )

    for sample, with_insertion in (("positive", True), ("negative", False)):
        writer = Writer(os.path.join(out_dir, f"{sample}.sam"), sample, rng)
        for target in TARGETS:
            contig = target[0]
            ref_hap = Haplotype(contig, reference[contig])
            if with_insertion and target[3] == "target_insertion":
                alt_hap = Haplotype(contig, alt, insertion_at=INS_POS + TSD_LEN)
                fragments(rng, ref_hap, target, writer)
                fragments(rng, alt_hap, target, writer, duplicates=3)
            else:
                fragments(rng, ref_hap, target, writer)
                fragments(rng, ref_hap, target, writer)
        probe_edge_artefact(rng, reference, writer)
        writer.close()


if __name__ == "__main__":
    main(sys.argv[1])
