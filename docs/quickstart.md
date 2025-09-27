# Quickstart

## Infectiousness: TOST and TOIT

```python
import numpy as np
from sarscov2_infectiousness import InfectiousnessParams, TOST, TOIT

params = InfectiousnessParams()

# TOST
tost = TOST(params=params, rng_seed=42)
x = np.linspace(-10, 10, 1001)
pdf_tost = tost.pdf(x)

# TOIT
toit = TOIT(params=params, rng_seed=42)
y = np.linspace(0, 30, 1001)
pdf_toit = toit.pdf(y)

# Samples and generation time
samples_tost = tost.rvs(10_000)
samples_toit = toit.rvs(10_000)
gen_times = toit.generation_time(10_000)
```

## TN93 distances

```python
# Requires: conda install -c bioconda tn93
# and: pip install -e ".[distance]"
from sarscov2_infectiousness.distance import tn93_pairs

pairs = tn93_pairs("aligned_sequences.fasta", threshold=0.02)  # seqid1, seqid2, distance
```

## Transmission linkage probabilities

```python
import numpy as np
from sarscov2_infectiousness.linkage import estimate_linkage_probability

gd = np.arange(0, 6, dtype=float)     # 0..5 SNPs
td = np.linspace(-7, 7, 15)           # -7..+7 days

p = estimate_linkage_probability(
    genetic_distance=gd,
    sampling_interval=td,
    intermediate_generations=(0, 1),
    no_intermediates=8,
    num_simulations=8000,
    rng_seed=123,
)
# p is a vector aligned with broadcasting rules if using grids
```

See infectiousness.md, distance.md, linkage.md for details and caveats.
