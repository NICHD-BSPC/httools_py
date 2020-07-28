This documentation uses [sphinx](http://www.sphinx-doc.org) to buid the documentation.

The built documentation can be found at
http://NICHD-BSPC.github.io/httools_py. If you want to build a local copy of the
documentation:

- create an environment
- activate it
- install sphinx into it
- run the Makefile in `docs`


That is:

```bash
# Create env
conda env create -p httools_py-docs \
  --file requirements-docs.yaml

# activate it
source activate httools_py-docs

# build the docs
cd docs
make html
```

The locally-built docs will be in `docs/_build/html/index.html`.
