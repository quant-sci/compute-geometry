# Installation

`compute-geometry` requires **Python 3.12+**.

## Using uv

```bash
uv add compute-geometry
```

## Using pip

```bash
pip install compute-geometry
```

Once installed, the package is importable as `cgeom`:

```python
import cgeom

print(cgeom.Point(1.0, 2.0))
```

## From source

To work on the library or build the documentation locally, clone the repository
and install it in editable mode with the documentation dependencies:

```bash
git clone https://github.com/kleyt0n/compute-geometry
cd compute-geometry

uv sync
uv pip install sphinx furo myst-parser sphinx-copybutton
```

Then build the HTML docs:

```bash
cd docs
make html
# open _build/html/index.html
```
