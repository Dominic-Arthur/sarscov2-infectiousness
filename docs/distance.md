# TN93 pairwise genetic distance

We wrap the VEG/HIV‑TRACE `tn93` CLI (C/C++), which is fast and handles IUPAC ambiguities.

Install the CLI:
```bash
conda install -c bioconda tn93
```

Install Python extras:
```bash
pip install -e ".[distance]"
```

Usage
```python
from sarscov2_infectiousness.distance import tn93_pairs
pairs = tn93_pairs("aligned_sequences.fasta", threshold=0.02)  # DataFrame: seqid1,seqid2,distance
```

Notes
- Input must be an aligned multi‑FASTA (same length). For SARS‑CoV‑2, use MAFFT/nextalign or similar.
- `tn93` is single‑threaded. To parallelize across cores:
  - Split sequences into blocks and run block pairs concurrently using Python threads or GNU parallel.
  - Use `-s other.fasta` for cross‑block comparisons; avoid duplicate (i,j)/(j,i).

Parallel outline (block schedule)
- Create blocks block_0000.fasta, block_0001.fasta, …
- For each i:
  - Run `tn93 block_i.fasta` (within-block)
  - For each j > i: `tn93 block_i.fasta -s block_j.fasta` (cross-block)
- Concatenate CSV outputs (skip duplicate headers).

Caveats
- All-pairs output is O(n^2); for n=10k, this is ~50M pairs. Prefer thresholds (e.g., 0.015–0.02).
- Store inputs/outputs on SSD; avoid too many tiny files by merging outputs as you go.

If you need a ready-made parallel scheduler, ask and we’ll add it here.
