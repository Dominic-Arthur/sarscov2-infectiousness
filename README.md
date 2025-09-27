# sarscov2-infectiousness

Mechanistic infectiousness models and related tools for SARS‑CoV‑2:

- infectiousness: Variable infectiousness (E/P/I) model following Hart et al. (2021), with TOST and TOIT distributions.
- distance: Fast TN93 pairwise genetic distance wrapper around the VEG/HIV-TRACE `tn93` CLI.
- linkage: Probabilistic transmission linkage from genetic and temporal distance using Monte Carlo and optional Numba acceleration.

Reference
- Hart WS, Maini PK, Thompson RN (2021). High infectiousness immediately before COVID‑19 symptom onset highlights the importance of continued contact tracing. eLife 10:e65534. https://doi.org/10.7554/eLife.65534

## Install

Base package (infectiousness model only):
```bash
pip install -e .
```

Optional extras:
- Linkage (Numba JIT):
```bash
pip install -e ".[linkage]"
```
- TN93 wrappers (needs external CLI + pandas):
```bash
conda install -c bioconda tn93
pip install -e ".[distance]"
```

Python >= 3.9 is recommended.

## Quickstart

Infectiousness (TOST/TOIT)
```python
import numpy as np
from sarscov2_infectiousness import InfectiousnessParams, TOST, TOIT

params = InfectiousnessParams()

# TOST: onset to transmission
tost = TOST(params=params, rng_seed=42)
x = np.linspace(-10, 10, 1001)
pdf_tost = tost.pdf(x)

# TOIT: start of presymptomatic to transmission
toit = TOIT(params=params, rng_seed=42)
y = np.linspace(0, 30, 1001)
pdf_toit = toit.pdf(y)

# Sampling
samples_tost = tost.rvs(10_000)
samples_toit = toit.rvs(10_000)
gen_times = toit.generation_time(10_000)
```

TN93 pairwise distances (thresholded)
```python
# Requires: conda install -c bioconda tn93; pip install ".[distance]"
from sarscov2_infectiousness.distance import tn93_pairs
pairs = tn93_pairs("aligned_sequences.fasta", threshold=0.02)  # DataFrame: seqid1, seqid2, distance
```

Transmission linkage probability
```python
import numpy as np
from sarscov2_infectiousness.linkage import estimate_linkage_probability

gd = np.array([0, 1, 2], dtype=float)          # SNP differences
td = np.array([-5.0, 0.0, 5.0], dtype=float)   # sampling interval (days)

p = estimate_linkage_probability(
    genetic_distance=gd,
    sampling_interval=td,
    intermediate_generations=(0, 1),  # include direct and single-intermediate links
    no_intermediates=10,              # simulate up to M=10 intermediates
    num_simulations=5000,
    rng_seed=1,
)
print(p)   # array of link probabilities in [0,1]
```

## Project structure

```
src/sarscov2_infectiousness/
  __init__.py
  infectiousness/
    __init__.py
    variable_model.py         # InfectiousnessParams, TOST, TOIT (+ presymptomatic_fraction helper)
  distance/
    __init__.py
    tn93.py                   # tn93_pairs wrapper
  linkage/
    __init__.py
    transmission_linkage_model.py
```

## Development

- Editable install:
  ```bash
  pip install -e ".[dev]"
  ```
- Tests:
  ```bash
  pytest -q
  ```
- Lint/format (optional):
  ```bash
  ruff check src tests
  black src tests
  ```

## Notes

- TOIT.pdf uses a vectorized trapezoidal integral; increase y_grid_points if you need higher accuracy.
- rvs() for TOIT/TOST uses a discretized grid; suitable for simulation but is an approximation.
- The `tn93` CLI is single‑threaded; for large datasets parallelize by running multiple instances across blocks (see docs/distance.md).

## License

MIT (see LICENSE)

## Citation

Please cite Hart et al. (2021) for the infectiousness model and acknowledge this package in analyses and publications.
