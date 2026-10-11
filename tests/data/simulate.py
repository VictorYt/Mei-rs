"""Simulate the synthetic test fixtures (see README.md).

Writes, into the directory given as the only argument:
- reference.fa: two random contigs; chr2 holds a reference copy of the Alu;
- targets.bed: five capture targets;
- targets_tiled.bed: the same, with the insertion target split into two
  overlapping probes whose edges fall on the insertion junctions;
- targets_tract.bed: one target around a reference poly(A) tract on chr2;
- positive.sam: heterozygous Alu, LINE-1 and SVA insertions with a TSD and a
  poly(A) tail on chr1, a probe-edge artefact, a few marked PCR duplicates,
  a heterozygous expansion of the chr2 poly(A) tract, and a chimeric
  junction with Alu clips and random unmapped mates on chr2;
- negative.sam: same targets, without the insertions;
- homozygous.sam: the Alu insertion on both haplotypes, in its target only;
- truth.tsv: the simulated insertions.

The elements are Dfam consensus sequences (data/consensus/mei_consensus.fa,
CC0), without their poly(A) tail. Reads are paired-end 2x100 bp and are
"aligned" from the simulation truth: anchored parts are matches, parts inside
an insertion are soft-clipped, reads entirely inside the Alu map to the chr2
copy with MAPQ 0, and the other reads inside an insertion are unmapped. Only
`random.random()` is used, so the output does not depend on the Python
version.
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

CONSENSUS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data", "consensus", "mei_consensus.fa")


def consensus(name):
    """A Dfam consensus sequence, without its poly(A) tail."""
    sequences, current = {}, None
    with open(CONSENSUS) as fasta:
        for line in fasta:
            line = line.strip()
            if line.startswith(">"):
                current = line[1:].split("#")[0]
                sequences[current] = []
            else:
                sequences[current].append(line)
    return "".join(sequences[name]).rstrip("A")


ALU = consensus("AluY")
POLYA = "A" * 30

CONTIGS = [("chr1", 20000), ("chr2", 10000)]


class Insertion:
    """An insertion after the 0-based position `pos` of chr1: the target site
    chr1[pos:pos + tsd_len] is duplicated, and the element (with its tail)
    is inserted in the given orientation:
    alt = chr1[:pos + tsd_len] + sequence + chr1[pos:], where `sequence` is
    the element + tail, reverse-complemented for the minus strand."""

    def __init__(self, name, element, pos, tsd_len, strand, target, copy=None):
        self.name = name
        self.element = element
        self.pos = pos
        self.tsd_len = tsd_len
        self.strand = strand
        forward = element + POLYA
        self.sequence = forward if strand == "+" else forward.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        self.target = target
        # Reference copy of the element (contig, position), where reads
        # entirely inside it map with MAPQ 0.
        self.copy = copy


# Reference Alu copy on chr2, outside every target.
ALU_COPY_POS = 7000
INSERTIONS = [
    Insertion("AluY", ALU, 10000, 12, "+", "target_insertion", ("chr2", ALU_COPY_POS)),
    # 5'-truncated LINE-1: the last 500 bp of the L1HS 3' end.
    Insertion("L1HS_3end", consensus("L1HS_3end")[-500:], 15000, 15, "-", "target_line1"),
    # The 3' end of an SVA_F.
    Insertion("SVA_F", consensus("SVA_F")[-600:], 6000, 10, "+", "target_sva"),
]
# The Alu insertion, used by the tiled targets below.
INS_POS = INSERTIONS[0].pos
TSD_LEN = INSERTIONS[0].tsd_len
# Chimeric junction of `positive` (random clips and mates), in target_negative.
CHIMERA_JUNCTION = 4150
# Probe whose start produces artefactual soft-clips.
PROBE_EDGE = ("chr1", 3000)
TARGETS = [
    ("chr1", 9850, 10150, "target_insertion"),
    ("chr1", 3000, 3200, "target_probe_edge"),
    ("chr2", 4000, 4300, "target_negative"),
    ("chr1", 14850, 15150, "target_line1"),
    ("chr1", 5850, 6150, "target_sva"),
]
# Same targets, the insertion one tiled with two overlapping probes: the first
# ends at the left junction, the second starts at the right one. Every
# soft-clip of the insertion stops at a probe edge, yet its poly(A) tail and
# TSD must keep it from the probe-edge filter.
TILED_TARGETS = [
    ("chr1", 9850, INS_POS + TSD_LEN, "target_insertion_left"),
    ("chr1", INS_POS, 10150, "target_insertion_right"),
] + TARGETS[1:]

# Reference poly(A) tract at the end of an Alu fragment on chr2 (outside the
# targets above), expanded by TRACT_EXPANSION bases on one haplotype of
# `positive`: alt = chr2[:TRACT_POS + TRACT_LEN] + A(expansion) + chr2[TRACT_POS + TRACT_LEN:].
TRACT_POS = 2000
TRACT_LEN = 12
TRACT_EXPANSION = 14
TRACT_FLANK = ALU[-40:]
TRACT_TARGET = ("chr2", 1850, 2150, "target_tract")

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

    def __init__(self, contig, reference, insertion=None):
        self.contig = contig
        self.insertion = insertion
        if insertion is None:
            self.sequence = reference
        else:
            ins = insertion
            self.sequence = reference[: ins.pos + ins.tsd_len] + ins.sequence + reference[ins.pos :]

    def to_alt(self, ref_pos):
        """Haplotype index of a reference position (second TSD copy side)."""
        ins = self.insertion
        if ins is None or ref_pos < ins.pos + ins.tsd_len:
            return ref_pos
        return ref_pos + ins.tsd_len + len(ins.sequence)

    def align(self, start, end):
        """Map the haplotype interval [start, end) to the reference."""
        ins = self.insertion
        if ins is None:
            return Mapping(self.contig, start, f"{end - start}M", 60)
        ins_start = ins.pos + ins.tsd_len
        ins_end = ins_start + len(ins.sequence)
        shift = len(ins.sequence) + ins.tsd_len
        if end <= ins_start:
            return Mapping(self.contig, start, f"{end - start}M", 60)
        if start >= ins_end:
            return Mapping(self.contig, start - shift, f"{end - start}M", 60)
        left = ins_start - start
        right = end - ins_end
        if left >= MIN_ANCHOR:
            return Mapping(self.contig, start, f"{left}M{end - start - left}S", 60)
        if right >= MIN_ANCHOR:
            return Mapping(self.contig, ins.pos, f"{end - start - right}S{right}M", 60)
        if ins.copy and start >= ins_start and end <= ins_start + len(ins.element):
            # Multi-mapping read: placed on the reference copy, MAPQ 0.
            contig, pos = ins.copy
            return Mapping(contig, pos + start - ins_start, f"{end - start}M", 0)
        return Mapping()


class TractHaplotype:
    """A haplotype of chr2 with the poly(A) tract expanded, aligned like a
    read aligner does: a read crossing the expansion is anchored on its longer
    side. Anchored upstream, it matches through the reference tract and the
    extra bases and the downstream flank are clipped (clip after the tract);
    anchored downstream, it matches back through the reference tract and the
    upstream flank and the extra bases are clipped (clip before the tract)."""

    def __init__(self, contig, reference):
        self.contig = contig
        self.end = TRACT_POS + TRACT_LEN
        self.sequence = reference[: self.end] + "A" * TRACT_EXPANSION + reference[self.end :]

    def to_alt(self, ref_pos):
        return ref_pos if ref_pos < self.end else ref_pos + TRACT_EXPANSION

    def align(self, start, end):
        shift = TRACT_EXPANSION
        if end <= self.end:
            return Mapping(self.contig, start, f"{end - start}M", 60)
        if start >= TRACT_POS + shift:
            # Starts in the tract: its A run matches the reference tract.
            return Mapping(self.contig, start - shift, f"{end - start}M", 60)
        upstream = self.end - start
        # Matches from the start of the tract once the extra bases are clipped.
        downstream = end - (TRACT_POS + shift)
        if upstream >= downstream and upstream >= MIN_ANCHOR:
            return Mapping(self.contig, start, f"{upstream}M{end - start - upstream}S", 60)
        if downstream >= MIN_ANCHOR:
            clip = end - start - downstream
            return Mapping(self.contig, TRACT_POS, f"{clip}S{downstream}M", 60)
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

    def pair(self, left, right, duplicate=False, name=None):
        """Write a fragment. `left` and `right` are (mapping, sequence) on the
        forward strand of the haplotype; `right` is sequenced in reverse."""
        (lmap, lseq), (rmap, rseq) = left, right
        if not lmap.mapped and not rmap.mapped:
            return
        rmap.reverse = True
        self.count += 1
        name = name or f"{self.sample}_{self.count:05d}"
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


def chimeras(rng, reference, writer, count=6):
    """A chimeric junction in `target_negative` (chr2:4150): reads soft-clipped
    with the first 25 bases of the Alu (typed), and fragments whose mate is
    unmapped and random (chimeric or adapter sequence). Its mates are checked
    (it is typed and has junction clips) and none is a mobile element."""
    seq = reference["chr2"]
    junction = CHIMERA_JUNCTION
    for i in range(3):
        start = junction - (READ_LEN - 25)
        anchored = seq[start:junction]
        left = (Mapping("chr2", start, f"{READ_LEN - 25}M25S", 60), anchored + ALU[:25])
        rstart = start + 150 + rng.below(50)
        right = (Mapping("chr2", rstart, f"{READ_LEN}M", 60), sequencing_errors(rng, seq[rstart : rstart + READ_LEN]))
        writer.pair(left, right, name=f"{writer.sample}_chimera_clip_{i + 1}")
    for i in range(count):
        start = junction - READ_LEN - rng.below(60)
        left = (Mapping("chr2", start, f"{READ_LEN}M", 60), sequencing_errors(rng, seq[start : start + READ_LEN]))
        right = (Mapping(), rng.sequence(READ_LEN))
        writer.pair(left, right, name=f"{writer.sample}_chimera_{i + 1}")


def main(out_dir):
    rng = Rng(SEED)
    reference = {name: rng.sequence(length) for name, length in CONTIGS}
    chr2 = reference["chr2"]
    chr2 = chr2[:ALU_COPY_POS] + ALU + chr2[ALU_COPY_POS + len(ALU) :]
    tract = TRACT_FLANK + "A" * TRACT_LEN
    start = TRACT_POS - len(TRACT_FLANK)
    reference["chr2"] = chr2[:start] + tract + chr2[start + len(tract) :]

    with open(os.path.join(out_dir, "reference.fa"), "w") as fa:
        for name, _ in CONTIGS:
            fa.write(f">{name}\n")
            seq = reference[name]
            for i in range(0, len(seq), 60):
                fa.write(seq[i : i + 60] + "\n")

    for file_name, targets in (
        ("targets.bed", TARGETS),
        ("targets_tiled.bed", TILED_TARGETS),
        ("targets_tract.bed", [TRACT_TARGET]),
    ):
        with open(os.path.join(out_dir, file_name), "w") as bed:
            for contig, start, end, name in targets:
                bed.write(f"{contig}\t{start}\t{end}\t{name}\n")

    chr1 = reference["chr1"]
    families = {"AluY": "ALU", "L1HS_3end": "LINE1", "SVA_F": "SVA"}
    with open(os.path.join(out_dir, "truth.tsv"), "w") as truth:
        truth.write(
            "#contig\tleft_junction\tright_junction\ttsd_start\ttsd_end\ttsd\telement\tfamily"
            "\tlength\tpolya\tstrand\tgenotype\n"
        )
        for ins in sorted(INSERTIONS, key=lambda i: i.pos):
            tsd = chr1[ins.pos : ins.pos + ins.tsd_len]
            truth.write(
                f"chr1\t{ins.pos + ins.tsd_len}\t{ins.pos + 1}\t{ins.pos + 1}\t{ins.pos + ins.tsd_len}\t{tsd}"
                f"\t{ins.name}\t{families[ins.name]}\t{len(ins.element)}\t{len(POLYA)}\t{ins.strand}\t0/1\n"
            )

    for sample, with_insertion in (("positive", True), ("negative", False)):
        writer = Writer(os.path.join(out_dir, f"{sample}.sam"), sample, rng)
        for target in TARGETS:
            contig = target[0]
            ref_hap = Haplotype(contig, reference[contig])
            insertion = next((i for i in INSERTIONS if i.target == target[3]), None)
            if with_insertion and insertion:
                alt_hap = Haplotype(contig, reference[contig], insertion)
                fragments(rng, ref_hap, target, writer)
                duplicates = 3 if insertion is INSERTIONS[0] else 0
                fragments(rng, alt_hap, target, writer, duplicates=duplicates)
            else:
                fragments(rng, ref_hap, target, writer)
                fragments(rng, ref_hap, target, writer)
        probe_edge_artefact(rng, reference, writer)
        if with_insertion:
            # Own random generator: the reads above do not change.
            main_rng, writer.rng = writer.rng, Rng(SEED + 1)
            contig = TRACT_TARGET[0]
            fragments(writer.rng, Haplotype(contig, reference[contig]), TRACT_TARGET, writer)
            fragments(writer.rng, TractHaplotype(contig, reference[contig]), TRACT_TARGET, writer)
            writer.rng = Rng(SEED + 2)
            chimeras(writer.rng, reference, writer)
            writer.rng = main_rng
        writer.close()

    # Homozygous Alu insertion (on both haplotypes) in target_insertion, with
    # its own random generator: the samples above do not change.
    rng = Rng(SEED + 3)
    writer = Writer(os.path.join(out_dir, "homozygous.sam"), "homozygous", rng)
    contig = TARGETS[0][0]
    alt_hap = Haplotype(contig, reference[contig], INSERTIONS[0])
    fragments(rng, alt_hap, TARGETS[0], writer)
    fragments(rng, alt_hap, TARGETS[0], writer)
    writer.close()


if __name__ == "__main__":
    main(sys.argv[1])
