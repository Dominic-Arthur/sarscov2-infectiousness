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
