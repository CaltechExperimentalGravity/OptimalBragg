SPHINXBUILD   = sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

.PHONY: help codemap html

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)"

# Render codemap from Graphviz source — keeps .dot and .svg/.png in sync
codemap: codemap.svg codemap.png

codemap.svg: codemap.dot
	dot -Tsvg $< -o $@

codemap.png: codemap.dot
	dot -Tpng $< -o $@ -Gdpi=150

# Build Sphinx HTML; render codemap first so images are fresh
html: codemap
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)"

# All other Sphinx targets (latexpdf, etc.)
%:
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)"
