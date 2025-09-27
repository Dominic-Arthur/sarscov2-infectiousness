#!/usr/bin/env python3
"""
Upgrade sarscov2-infectiousness to a subpackage layout:

src/sarscov2_infectiousness/
  infectiousness/variable_model.py   (moved from top-level)
  distance/tn93.py                   (new)
  linkage/transmission_linkage_model.py (new)
  __init__.py                        (updated)

Also updates pyproject.toml to add optional extras:
  [project.optional-dependencies]
  linkage = ["numba>=0.59"]
  distance = ["pandas>=2.0"]

Run this script from the repository root.
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

ROOT = Path(".").resolve()
SRC_ROOT = ROOT / "src" / "sarscov2_infectiousness"

INF_DIR = SRC_ROOT / "infectiousness"
DIST_DIR = SRC_ROOT / "distance"
LINK_DIR = SRC_ROOT / "linkage"

OLD_VM = SRC_ROOT / "variable_model.py"
NEW_VM = INF_DIR / "variable_model.py"

TOP_INIT = SRC_ROOT / "__init__.py"
INF_INIT = INF_DIR / "__init__.py"
DIST_INIT = DIST_DIR / "__init__.py"
LINK_INIT = LINK_DIR / "__init__.py"

PYPROJECT = ROOT / "pyproject.toml"

TESTS_DIR = ROOT / "tests"
TEST_LINKAGE = TESTS_DIR / "test_linkage.py"


def ensure_dirs():
    for d in [INF_DIR, DIST_DIR, LINK_DIR, TESTS_DIR]:
        d.mkdir(parents=True, exist_ok=True)


def move_variable_model():
    if OLD_VM.exists():
        # Back up old file before moving
        backup = OLD_VM.with_suffix(".py.bak")
        if not backup.exists():
            shutil.copy2(OLD_VM, backup)
            print(f"Backed up {OLD_VM} -> {backup}")
        # Move into infectiousness subpackage
        shutil.move(str(OLD_VM), str(NEW_VM))
        print(f"Moved {OLD_VM} -> {NEW_VM}")
    elif not NEW_VM.exists():
        # If no existing file, write a minimal placeholder (should not happen if you used the scaffold)
        NEW_VM.write_text(
            "# Placeholder: variable_model.py not found; please paste your model here.\n",
            encoding="utf-8",
        )
        print(f"Wrote placeholder {NEW_VM}")
    else:
        print(f"variable_model already in place at {NEW_VM}")


def write_infectiousness_init():
    content = '''\
from .variable_model import InfectiousnessParams, TOST, TOIT
# presymptomatic_fraction may exist in your variable_model; import if present.
try:
    from .variable_model import presymptomatic_fraction
except Exception:  # pragma: no cover
    presymptomatic_fraction = None  # type: ignore

__all__ = ["InfectiousnessParams", "TOST", "TOIT"]
if presymptomatic_fraction is not None:
    __all__.append("presymptomatic_fraction")
'''
    INF_INIT.write_text(content, encoding="utf-8")
    print(f"Wrote {INF_INIT}")


def write_top_init():
    content = '''\
"""
sarscov2_infectiousness: infectiousness models and related tools.
Keep top-level imports light to avoid hard optional dependencies.
"""

from .infectiousness import InfectiousnessParams, TOST, TOIT

# Export helper if present
try:
    from .infectiousness import presymptomatic_fraction  # type: ignore
except Exception:  # pragma: no cover
    presymptomatic_fraction = None  # type: ignore

__all__ = ["InfectiousnessParams", "TOST", "TOIT"]
if presymptomatic_fraction is not None:
    __all__.append("presymptomatic_fraction")

__version__ = "0.1.0"
'''
    # Backup existing init if different
    if TOP_INIT.exists():
        orig = TOP_INIT.read_text(encoding="utf-8")
        if orig != content:
            backup = TOP_INIT.with_suffix(".py.bak")
            if not backup.exists():
                shutil.copy2(TOP_INIT, backup)
                print(f"Backed up {TOP_INIT} -> {backup}")
    TOP_INIT.write_text(content, encoding="utf-8")
    print(f"Wrote {TOP_INIT}")


def write_distance():
    DIST_INIT.write_text('from .tn93 import tn93_pairs\n__all__ = ["tn93_pairs"]\n', encoding="utf-8")
    tn93_py = DIST_DIR / "tn93.py"
    tn93_py.write_text(
        '''\
from __future__ import annotations
import io
import shutil
import subprocess
from pathlib import Path
import pandas as pd

def _ensure_tn93_available() -> None:
    if shutil.which("tn93") is None:
        raise RuntimeError(
            "The 'tn93' CLI was not found on PATH. Install with conda: "
            "conda install -c bioconda tn93"
        )

def tn93_pairs(fasta_path: str | Path, threshold: float = 0.02, quiet: bool = True) -> pd.DataFrame:
    """
    Run VEG/HIV-TRACE tn93 on a multi-FASTA alignment and return a pandas DataFrame:
    columns: seqid1, seqid2, distance. Threshold filters pairs (default 0.02).
    """
    _ensure_tn93_available()
    cmd = ["tn93", "-t", str(threshold), "-f", "csv", str(fasta_path)]
    if quiet:
        cmd.append("-q")
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return pd.read_csv(io.StringIO(result.stdout))
''',
        encoding="utf-8",
    )
    print(f"Wrote {DIST_INIT} and {tn93_py}")


def write_linkage():
    LINK_INIT.write_text(
        'from .transmission_linkage_model import estimate_linkage_probability, pairwise_linkage_probability_matrix\n'
        '__all__ = ["estimate_linkage_probability", "pairwise_linkage_probability_matrix"]\n',
        encoding="utf-8",
    )
    tlm_py = LINK_DIR / "transmission_linkage_model.py"
    tlm_py.write_text(
        '''\
"""
Estimate the probability of a transmission link between two cases from
genetic distance (SNPs) and temporal distance (days).

Depends on:
- sarscov2_infectiousness.infectiousness (TOIT/params)
- Optional: numba for JIT acceleration; falls back to pure-Python if unavailable.

Public API:
- estimate_linkage_probability(...)
- pairwise_linkage_probability_matrix(...)
"""

from __future__ import annotations

from typing import Iterable, Tuple, Union, Dict

import numpy as np

# Import infectiousness model
from sarscov2_infectiousness.infectiousness import TOIT, InfectiousnessParams

ArrayLike = Union[float, Iterable[float], np.ndarray]

# Optional Numba: if not installed, provide a no-op decorator
try:
    import numba  # type: ignore
    def njit(fn=None, **kw):
        if fn is None:
            return lambda f: numba.njit(cache=True,(
        tmrca: np.ndarray,  # shape (N,)
        tmrca_expected: np.ndarray,  # shape (N,)
        tolerance: np.ndarray  # shape (N,)
) -> float:
    """
    Mean over simulations of the event:
    |tmrca_j - tmrca_expected_j| <= tolerance_j
    """
    N = tmrca.shape[0]
    count = 0
    for j in range(N):
        if abs(tmrca[j] - tmrca_expected[j]) <= tolerance[j]:
            count += 1
    return count / N

@njit
def _genetic_kernel_many(
        dists: np.ndarray,  # shape (K,)
        clock_rates: np.ndarray,  # shape (N,)
        psi_sA: np.ndarray,  # shape (N,)
        psi_sB: np.ndarray,  # shape (N,)
        generation_interval: np.ndarray,  # shape (N, M+1)
        toit_difference: np.ndarray,  # shape (N,)
        incubation_period: np.ndarray,  # shape (N, 2)
        caseX_to_caseA: np.ndarray,  # shape (N,)
        intermediates: int  # M (max number of intermediates)
) -> np.ndarray:
    """
    For each genetic distance d in dists, compute vector p_m for m=0..M,
    where p_m is the genetic evidence under the scenario of m intermediates.
    Returns an array of shape (K, M+1).
    """
    K = dists.shape[0]
    N = clock_rates.shape[0]
    M = intermediates

    out = np.zeros((K, M + 1), dtype=np.float64)

    # Invariants across m
    dir_tmrca_exp = psi_sA + psi_sB  # shape (N,)
    inc_period_sum = incubation_period[:, 0] + incubation_period[:, 1]  # shape (N,)

    for k in range(K):
        d = dists[k]
        # tmrca per simulation for this distance
        tmrca = d / (2.0 * clock_rates)  # shape (N,)

        # m = 0: direct
        p_direct = _mean_prob_kernel(tmrca, dir_tmrca_exp, generation_interval[:, 0])
        out[k, 0] = p_direct

        # m > 0
        # The original indexing logic is preserved for consistency
        for m in range(1, M + 1):
            idx = M - (m - 1)

            # Successive transmission path
            # sum generation intervals from idx to end (inclusive)
            suc_sum = np.zeros(N, dtype=np.float64)
            for j in range(N):
                s = 0.0
                for col in range(idx, M + 1):
                    s += generation_interval[j, col]
                suc_sum[j] = s
            suc_tmrca = dir_tmrca_exp + suc_sum
            p_suc = _mean_prob_kernel(tmrca, suc_tmrca, generation_interval[:, 0])

            # Common ancestor path: sum from idx to M-1 (exclusive of last col)
            com_sum = np.zeros(N, dtype=np.float64)
            for j in range(N):
                s = 0.0
                for col in range(idx, M):
                    s += generation_interval[j, col]
                com_sum[j] = s
            common_tmrca = toit_difference + inc_period_sum + com_sum
            p_common = _mean_prob_kernel(tmrca, common_tmrca, caseX_to_caseA)

            # Inclusion-exclusion
            out[k, m] = p_common + p_suc - (p_common * p_suc)

    return out


def _run_simulations(toit: TOIT, num_simulations: int, no_intermediates: int) -> Dict[str, np.ndarray]:
    """
    Draw all random epidemiological quantities once.
    """
    N = int(num_simulations)
    M = int(no_intermediates)

    inc_period = toit.sample_incubation(size=(N, 2))  # (N, 2)
    gen_interval = toit.generation_time(size=(N, M + 1))  # (N, M+1)
    toit_values = toit.rvs(size=(N, 2))  # (N, 2)
    latent_period = toit.sample_E(size=N)  # (N,)
    clock_rates = toit.sample_clock_rate_per_day(size=N)  # (N,)

    sim = {
        "incubation_period": inc_period,
        "generation_interval": gen_interval,
        "psi_sA": np.abs(gen_interval[:, 0] - inc_period[:, 0]),
        "psi_sB": inc_period[:, 1],
        "diff_inc": inc_period[:, 0] - inc_period[:, 1],
        "caseX_to_caseA": latent_period + np.minimum(toit_values[:, 0], toit_values[:, 1]),
        "toit_difference": np.abs(toit_values[:, 1] - toit_values[:, 0]),
        "clock_rates": clock_rates,
    }
    return sim


def estimate_linkage_probability(
        genetic_distance: ArrayLike,  # SNPs; scalar or array
        sampling_interval: ArrayLike,  # days; scalar or array (same length as genetic_distance)
        intermediate_generations: Tuple[int, ...] = (0,),  # which m to include in final mixture
        no_intermediates: int = 10,  # max intermediates M used in simulation/kernel
        infectiousness_profile: InfectiousnessParams | None = None,  # Use default InfectiousnessParams
        subs_rate: float = 1e-3,  # subs/site/year (median)
        subs_rate_sigma: float = 0.33,  # lognormal sigma for relaxed clock
        relax_rate: bool = False,  # relaxed clock on/off
        num_simulations: int = 10000,  # Monte Carlo draws
        rng_seed: int = 12345,  # passed into TOIT
) -> Union[float, np.ndarray]:
    """
    Estimate P(link | g, t) combining temporal and genetic evidence.
    Returns float if both inputs are scalars; else a 1D array.
    """
    # 1) Prepare inputs as 1D arrays of same length
    g = np.atleast_1d(np.asarray(genetic_distance, dtype=float))
    t = np.atleast_1d(np.asarray(sampling_interval, dtype=float))
    if g.shape[0] != t.shape[0]:
        raise ValueError(
            f"genetic_distance and sampling_interval must have the same length, got {g.shape[0]} vs {t.shape[0]}."
        )

    K = g.shape[0]
    M = int(no_intermediates)

    # 2) Initialize model and run simulations once
    toit = TOIT(
        params=infectiousness_profile,
        rng_seed=int(rng_seed),
        subs_rate=float(subs_rate),
        subs_rate_sigma=float(subs_rate_sigma),
        relax_rate=bool(relax_rate),
    )
    sim = _run_simulations(toit, int(num_simulations), M)

    # 3) Temporal evidence
    p_t = _temporal_kernel(
        sampling_interval=t,
        diff_inc=sim["diff_inc"],
        generation_interval0=sim["generation_interval"][:, 0],
    )  # shape (K,)

    # 4) Genetic evidence: p_m for each m=0..M
    p_m = _genetic_kernel_many(
        dists=g,
        clock_rates=sim["clock_rates"],
        psi_sA=sim["psi_sA"],
        psi_sB=sim["psi_sB"],
        generation_interval=sim["generation_interval"],
        toit_difference=sim["toit_difference"],
        incubation_period=sim["incubation_period"],
        caseX_to_caseA=sim["caseX_to_caseA"],
        intermediates=M,
    )  # shape (K, M+1)

    # 5) Normalize genetic evidence across m
    total = 1.0 - np.prod(1.0 - p_m, axis=1)  # shape (K,)
    with np.errstate(divide="ignore", invalid="ignore"):
        p_rel = np.where(total[:, None] > 0.0, p_m / total[:, None], 0.0)
    row_sums = p_rel.sum(axis=1)  # shape (K,)
    with np.errstate(divide="ignore", invalid="ignore"):
        p_norm = np.where(row_sums[:, None] > 0.0, p_rel / row_sums[:, None], 0.0)

    # Select columns m specified by intermediate_generations
    cols = np.array(intermediate_generations, dtype=np.int64)
    if cols.min() < 0 or cols.max() > M:
        raise ValueError(f"intermediate_generations must be within [0, {M}], got {intermediate_generations}.")
    selected = p_norm[:, cols].sum(axis=1)  # shape (K,)

    # 6) Combine temporal and genetic evidence
    out = p_t * selected  # shape (K,)

    # 7) Return scalar if scalar input
    if np.isscalar(genetic_distance) and np.isscalar(sampling_interval):
        return float(out[0])
    return out


def pairwise_linkage_probability_matrix(
        genetic_distances: np.ndarray,  # 1D
        temporal_distances: np.ndarray,  # 1D
        intermediate_generations: Tuple[int, ...] = (0,),
        no_intermediates: int = 10,
        **kwargs,
) -> np.ndarray:
    """
    Compute a (len(genetic_distances) x len(temporal_distances)) matrix of P(link | g, t).
    """
    gd = np.asarray(genetic_distances, dtype=float)
    td = np.asarray(temporal_distances, dtype=float)

    g_grid, t_grid = np.meshgrid(gd, td, indexing="ij")  # each shape (G, T)
    flat_p = estimate_linkage_probability(
        genetic_distance=g_grid.ravel(),
        sampling_interval=t_grid.ravel(),
        intermediate_generations=intermediate_generations,
        no_intermediates=no_intermediates,
        **kwargs,
    )
    return np.asarray(flat_p, dtype=float).reshape(g_grid.shape)
''',
        encoding="utf-8",
    )
    print(f"Wrote {LINK_INIT} and {tlm_py}")


def update_pyproject():
    if not PYPROJECT.exists():
        print("pyproject.toml not found; skipping update.")
        return

    text = PYPROJECT.read_text(encoding="utf-8")

    # Ensure optional-dependencies section exists
    if "[project.optional-dependencies]" not in text:
        text += "\n[project.optional-dependencies]\n"

    # Add/ensure linkage extra
    if re.search(r"^\s*linkage\s*=", text, flags=re.MULTILINE) is None:
        text = re.sub(
            r"(\[project\.optional-dependencies\]\s*)",
            r"\\1linkage = [\"numba>=0.59\"]\n",
            text,
            count=1,
            flags=re.MULTILINE,
        )

    # Add/ensure distance extra
    if re.search(r"^\s*distance\s*=", text, flags=re.MULTILINE) is None:
        text = re.sub(
            r"(\[project\.optional-dependencies\]\s*(?:.*\n)*)$",
            r"\1distance = [\"pandas>=2.0\"]\n",
            text,
            count=1,
            flags=re.DOTALL,
        )

    PYPROJECT.write_text(text, encoding="utf-8")
    print("Updated pyproject.toml (added optional-dependencies: linkage, distance)")


def write_linkage_test():
    TESTS_DIR.mkdir(parents=True, exist_ok=True)
    TEST_LINKAGE.write_text(
        '''\
import numpy as np
from sarscov2_infectiousness.linkage import estimate_linkage_probability

def test_linkage_shapes():
    gd = np.array([0, 1, 2], dtype=float)
    td = np.array([-5.0, 0.0, 5.0], dtype=float)
    p = estimate_linkage_probability(
        genetic_distance=gd,
        sampling_interval=td,
        num_simulations=500,  # keep test fast
        no_intermediates=3,
        intermediate_generations=(0,1),
        rng_seed=1,
    )
    assert p.shape == (3,)
    assert np.all((p >= 0) & (p <= 1))
''',
        encoding="utf-8",
    )
    print(f"Wrote {TEST_LINKAGE}")


def main():
    # Sanity checks
    if not (ROOT / "pyproject.toml").exists():
        print("Please run this script from the repository root (pyproject.toml not found).", file=sys.stderr)
        sys.exit(1)
    if not SRC_ROOT.exists():
        print("src/sarscov2_infectiousness not found. Are you in the right project?", file=sys.stderr)
        sys.exit(1)

    ensure_dirs()
    move_variable_model()
    write_infectiousness_init()
    write_top_init()
    write_distance()
    write_linkage()
    update_pyproject()
    write_linkage_test()

    print("\nUpgrade complete.")
    print("Next steps:")
    print("  1) pip install -e '.[dev]'  # reinstall editable package")
    print("  2) pytest -q                # run tests")
    print("  3) If you plan to use tn93 wrappers, install tn93 via conda and pandas via:")
    print("       conda install -c bioconda tn93")
    print("       pip install '.[distance]'")
    print("  4) For linkage acceleration, install numba via:")
    print("       pip install '.[linkage]'")

if __name__ == "__main__":
    main()
