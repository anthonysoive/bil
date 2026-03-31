import os
import sys

project = "Bil"
copyright = "Université Gustave Eiffel"
author = "Patrick Dangla et al."
release = "2.0"

extensions = [
    "breathe",
    "myst_parser",
    "sphinx_copybutton",
    "sphinx.ext.mathjax",
    "sphinx.ext.autosectionlabel",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "furo"
html_static_path = ["_static"]
html_title = "Bil Documentation"

# Breathe config
# The Doxygen XML is generated to doc/oxygen/xml/ relative to the project root.
# conf.py lives in doc/sphinx/, so we go up two levels.
_doxy_xml = os.path.join(os.path.dirname(__file__), "..", "oxygen", "xml")
breathe_projects = {"Bil": _doxy_xml}
breathe_default_project = "Bil"
breathe_default_members = ("members", "undoc-members")

# MyST (Markdown) config
myst_enable_extensions = ["dollarmath", "colon_fence"]

autosectionlabel_prefix_document = True
