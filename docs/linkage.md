# Transmission linkage model

Estimate P(link | g, t) between two cases from:
- Genetic distance g (SNPs between consensus sequences)
- Temporal distance t (days between sample dates)

Approach
- Draw Monte Carlo samples of epidemiological quantities from the infectiousness model (TOIT, incubation).
- Model a relaxed or strict molecular clock to convert g to a time to MRCA (TMRCA).
- Compute two genetic evidence components per number of intermediates m (0..M):
  - Successive-transmission path
  - Common-ancestor path
- Combine via inclusion–exclusion and normalize across m.
- Multiply by temporal evidence P_t derived from Monte Carlo sampling of incubation/generation times.

API
- estimate_linkage_probability(genetic_distance, sampling_interval, intermediate_generations=(0,), no_intermediates=10, infectiousness_profile=None, subs_rate=1e-3, subs_rate_sigma=0.33, relax_rate=False, num_simulations=10000, rng_seed=12345)
  - Inputs:
    - genetic_distance: float or array (SNPs)
    - sampling_interval: float or array (days)
    - intermediate_generations: tuple of m to include (e.g., (0,) direct only; (0,1) direct or one intermediate)
    - no_intermediates: M, max m used in simulation/kernel
    - infectiousness_profile: InfectiousnessParams or None (defaults)
    - Clock parameters: subs/site/year (median), lognormal sigma, relaxed True/False
    - num_simulations: Monte Carlo draws (trade speed vs. variance)
    - rng_seed: reproducibility
  - Returns float if both inputs are scalars; else a vector aligned with broadcasting.

- pairwise_linkage_probability_matrix(genetic_distances, temporal_distances, …) → matrix

Example
```python
import numpy as np
from sarscov2_infectiousness.linkage import estimate_linkage_probability

gd = np.array([0, 1, 2, 3], dtype=float)
td = np.array([-5, 0, 5], dtype=float)

p = estimate_linkage_probability(
    genetic_distance=gd,
    sampling_interval=td,
    intermediate_generations=(0, 1),
    no_intermediates=8,
    relax_rate=True,         # use relaxed lognormal clock
    subs_rate=1e-3,
    subs_rate_sigma=0.33,
    num_simulations=8000,
    rng_seed=2,
)
```

Performance
- Install Numba (pip install ".[linkage]") for JIT acceleration of kernels.
- Complexity scales with:
  - num_simulations N (linear)
  - no_intermediates M (linear inside kernels)
  - number of query pairs K (linear)
- For large K, vectorize inputs where possible and reuse the same simulation outputs when evaluating many pairs under identical simulation settings.

Caveats
- Uses consensus SNP distances; within-host diversity is not explicitly modeled.
- Assumes an aligned genome and a simple molecular clock; choose subs_rate and sigma appropriate to your dataset.
- Outputs are probabilities under the model, not definitive transmission confirmation.

Reproducibility
- Set rng_seed for deterministic simulations.
- Fix infectiousness parameters and clock settings when comparing scenarios.
