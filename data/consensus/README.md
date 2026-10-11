# MEI consensus sequences

`mei_consensus.fa` holds the consensus sequences that `mei-rs scan` uses by
default to type the inserted element (built into the binary; `--consensus`
replaces them):

| Family | Sequences |
|---|---|
| Alu | 42 subfamilies: AluY (21, including AluYa5, AluYb8, AluYb9), AluS (17), AluJ (4) |
| LINE-1 | L1HS 5' end (2,136 bp) and 3' end (902 bp); L1PA2 to L1PA8 3' ends; L1P1 ORF2 (3,294 bp) |
| SVA | SVA_A to SVA_F |

They come from [Dfam](https://dfam.org) 4.0 (May 2026), released under the
[CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) license (public
domain dedication), and are downloaded by `fetch_dfam.py` (Python standard
library only):

```bash
python3 -I data/consensus/fetch_dfam.py data/consensus/mei_consensus.fa
```

Headers are `>name#class accession`, with the class `SINE/Alu`, `LINE/L1` or
`Retroposon/SVA` (RepeatMasker style); a custom `--consensus` file uses the
same form, or names starting with `Alu`, `L1` or `SVA`.

Dfam models L1HS and L1PA as their 5' and 3' ends only: the middle of the
element (most of ORF2, about 3 kb) comes from the ORF2 model of L1P1, the
youngest primate LINE-1 with one in Dfam. It is 96-100% identical to the
clips of young LINE-1 insertions whose 5' junction falls in ORF2 (#56);
without it they stay `untyped`. Such insertions are reported with the
subfamily `L1P1_orf2` and coordinates in that model.
