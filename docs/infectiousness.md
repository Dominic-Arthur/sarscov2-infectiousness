# Infectiousness model (E/P/I)

Implements the variable infectiousness model in Hart et al. (2021) for SARS‑CoV‑2, with Gamma-distributed stage durations:
- E (latent, not infectious)
- P (presymptomatic infectious)
- I (symptomatic infectious)

Notation (paper → code):
- k_inc, scale_inc: incubation Gamma shape/scale (E+P).
- k_E: latent shape; k_P = k_inc − k_E (derived).
- mu: symptomatic rate μ (mean symptomatic duration ≈ 1/μ when k_I=1).
- k_I: symptomatic shape (1 → exponential).
- alpha: relative infectiousness level P vs I.
- gamma_rate g = 1/(k_inc * scale_inc) (derived).
- C = (k_inc g μ) / (α k_P μ + k_inc g) (derived).

Distributions:
- TOST x = time from symptom onset to transmission (piecewise pdf):
  - x < 0: f(x) = α C [1 − F_P(−x)]
  - x ≥ 0: f(x) = C [1 − F_I(x)]
- TOIT y⋆ = time from start of P to transmission:
  - y⋆ ≥ 0: f(y⋆) = C [ α(1 − F_P(y⋆)) + ∫₀^{y⋆} (1 − F_I(y⋆ − y_P)) f_P(y_P) dy_P ]
  - Else 0.

API (summary)
- InfectiousnessParams: parameter container with derived properties k_P, scale_I, gamma_rate, C.
- TOST: pdf(x), rvs(size).
- TOIT: pdf(x), rvs(size), generation_time(size), sample_clock_rate_per_day(size).
- Helper: presymptomatic_fraction(params) returns q_P = (α k_P μ) / (α k_P μ + k_inc g).

Usage
```python
from sarscov2_infectiousness import InfectiousnessParams, TOIT, TOST

p = InfectiousnessParams()
toit = TOIT(params=p, rng_seed=1)
tost = TOST(params=p)

# Evaluate PDFs
import numpy as np
x = np.linspace(-10, 10, 5001); px = tost.pdf(x)
y = np.linspace(0, 30, 5001); py = toit.pdf(y)

# Numerical check (area ~ 1)
ax = np.trapz(px, x); ay = np.trapz(py, y)
```

Numerical notes
- TOIT.pdf uses a vectorized trapezoidal integral over y_P; increase y_grid_points for accuracy.
- rvs() uses discretized inverse transform on [a, b]; it’s an approximation to the continuous distribution.
- For generation_time, TOIT draws E + TOIT; treat as a proxy under the model assumptions.

Reference: Hart WS et al. (2021), eLife 10:e65534.
