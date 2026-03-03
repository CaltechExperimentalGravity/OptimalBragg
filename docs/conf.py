# Sphinx configuration for OptimalBragg documentation

project = 'OptimalBragg'
author = 'Caltech Experimental Gravity Group'
release = '0.2.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
]

# Napoleon settings (NumPy-style docstrings)
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True

# Autodoc settings
autodoc_member_order = 'bysource'
autodoc_default_options = {
    'members': True,
    'undoc-members': False,
    'show-inheritance': True,
}

# Theme — pydata-sphinx-theme supports light/dark mode switching
html_theme = 'pydata_sphinx_theme'
html_theme_options = {
    'navigation_depth': 3,
    'show_toc_level': 2,
}

# GitHub context for pydata theme
html_context = {
    'github_user': 'CaltechExperimentalGravity',
    'github_repo': 'OptimalBragg',
    'github_version': 'master',
    'doc_path': 'docs/',
}

# Exclude patterns
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
