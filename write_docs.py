#!/usr/bin/env python3
"""
Write or update README.md and docs/*.md for sarscov2-infectiousness.

Usage:
  - Place this script at the repository root (same level as pyproject.toml).
  - Run: python write_docs.py

Behavior:
  - Creates docs/ if missing.
  - Writes/updates:
      README.md
      docs/index.md
      docs/installation.md
      docs/quickstart.md
      docs/infectiousness.md
      docs/distance.md
      docs/linkage.md
      docs/api.md
  - If a file exists and content differs, a .bak backup is written once.
"""

from __future__ import annotations

import sys
from pathlib import Path
import textwrap
from typing import Dict


def write_file(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    new = content if content.endswith("\n") else content + "\n"
    if path.exists():
        old = path.read_text(encoding="utf-8")
        if old == new:
            print(f"Unchanged: {path}")
            return
        # Backup once
        bak = path.with_suffix(path.suffix + ".bak")
        if not bak.exists():
            bak.write_text(old, encoding="utf-8")
            print(f"Backed up: {path} -> {bak}")
    path.write_text(new, encoding="utf-8", newline="\n")
    print(f"Wrote:     {path}")


def main() -> None:
    root = Path(".").resolve()
    docs = root / "docs"

    files: Dict[Path, str] = {
        root / "README.md": textwrap.dedent(
            """\
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
            """
        ),
        docs / "index.md": textwrap.dedent(
            """\
            # sarscov2-infectiousness

            Mechanistic infectiousness models and tools for SARS‑CoV‑2:

            - Infectiousness: E/P/I stage model (TOST, TOIT) following Hart et al. (2021).
            - Distance: TN93 pairwise genetic distance via the VEG/HIV‑TRACE `tn93` CLI.
            - Linkage: Probability of transmission linkage from genetic and temporal distances.

            Get started:
            - Installation: installation.md
            - Quickstart: quickstart.md
            - Infectiousness model details: infectiousness.md
            - TN93 usage and scaling: distance.md
            - Linkage model: linkage.md
            - API: api.md
            """
        ),
        docs / "installation.md": textwrap.dedent(
            """\
            # Installation

            Base package (infectiousness only):
            ```bash
            pip install -e .
            ```

            Optional extras:
            - Linkage (Numba JIT acceleration):
            ```bash
            pip install -e ".[linkage]"
            ```
            - Distance (TN93 wrapper):
              - Install the external CLI:
                ```bash
                conda install -c bioconda tn93
                ```
              - Install Python extras:
                ```bash
                pip install -e ".[distance]"
                ```

            Python >= 3.9, NumPy >= 1.23, SciPy >= 1.9.

            Tip: Use a virtual environment (conda or venv). For Apple Silicon, prefer conda-forge builds for speed.
            """
        ),
        docs / "quickstart.md": textwrap.dedent(
            """\
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
            """
        ),
        docs / "infectiousness.md": textwrap.dedent(
            """\
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
            """
        ),
        docs / "distance.md": textwrap.dedent(
            """\
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
            """
        ),
        docs / "linkage.md": textwrap.dedent(
            """\
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
            """
        ),
        docs / "api.md": textwrap.dedent(
            """\
            # API

            Public modules and symbols:

            - sarscov2_infectiousness.infectiousness
              - InfectiousnessParams
              - TOST
              - TOIT
              - presymptomatic_fraction(params)  # if present

            - sarscov2_infectiousness.distance
              - tn93_pairs(fasta_path, threshold=0.02, quiet=True) → pandas.DataFrame

            - sarscov2_infectiousness.linkage
              - estimate_linkage_probability(...)
              - pairwise_linkage_probability_matrix(...)

            If you use MkDocs with mkdocstrings, you can auto-generate API docs from docstrings:

            1) Install:
            ```bash
            pip install mkdocs mkdocstrings[python] mkdocs-material
            ```

            2) mkdocs.yml:
            ```yaml
            site_name: sarscov2-infectiousness
            theme:
              name: material
            plugins:
              - mkdocstrings:
                  handlers:
                    python:
                      options:
                        docstring_style: numpy
            nav:
              - Home: index.md
              - Installation: installation.md
              - Quickstart: quickstart.md
              - Infectiousness: infectiousness.md
              - TN93 Distance: distance.md
              - Linkage: linkage.md
              - API: api.md
            ```

            3) docs/api.md contents (mkdocstrings blocks):
            ```md
            # API Reference

            ## Infectiousness
            ::: sarscov2_infectiousness.infectiousness.variable_model

            ## Distance
            ::: sarscov2_infectiousness.distance.tn93

            ## Linkage
            ::: sarscov2_infectiousness.linkage.transmission_linkage_model
            ```

            4) Serve:
            ```bash
            mkdocs serve
            ```
            """
        ),
    }

    for path, content in files.items():
        write_file(path, content)

    print("\nDone. Next steps:")
    print("  - Review README.md and docs/*.md")
    print("  - Optionally add mkdocs.yml and serve docs with `mkdocs serve`")
    print("  - Commit changes to your repository")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
