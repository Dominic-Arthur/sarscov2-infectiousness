# sarscov2-infectiousness

Variable infectiousness (E/P/I) model for SARS‑CoV‑2, implementing the TOST and TOIT distributions from Hart et al. (2021) [1].

## Install

```bash
pip install -e ".[dev]"
# or
pip install -e .
```

## Quickstart

```python
import numpy as np
from sarscov2_infectiousness import InfectiousnessParams, TOST, TOIT

params = InfectiousnessParams()  # defaults match [1] incubation; tune as needed

tost = TOST(params=params, rng_seed=42)
x = np.linspace(-10, 10, 1001)
pdf_tost = tost.pdf(x)

toit = TOIT(params=params, rng_seed=42)
y = np.linspace(0, 30, 1001)
pdf_toit = toit.pdf(y)

samples_tost = tost.rvs(10_000)
samples_toit = toit.rvs(10_000)
gen_times = toit.generation_time(10_000)
```

## Reference

- [1] Hart WS, Maini PK, Thompson RN (2021). High infectiousness immediately before COVID‑19 symptom onset highlights the importance of continued contact tracing. eLife 10:e65534. https://doi.org/10.7554/eLife.65534

## License

MIT (see LICENSE)
