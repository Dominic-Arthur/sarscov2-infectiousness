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
