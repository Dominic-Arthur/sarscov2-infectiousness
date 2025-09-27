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
