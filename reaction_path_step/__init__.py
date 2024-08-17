# -*- coding: utf-8 -*-

"""
reaction_path_step
A SEAMM plugin for finding transition states and reaction paths
"""

# Bring up the classes so that they appear to be directly in
# the reaction_path_step package.

from .reaction_path import ReactionPath  # noqa: F401, E501
from .reaction_path_parameters import ReactionPathParameters  # noqa: F401, E501
from .reaction_path_step import ReactionPathStep  # noqa: F401, E501
from .tk_reaction_path import TkReactionPath  # noqa: F401, E501

from .metadata import metadata  # noqa: F401

# Handle versioneer
from ._version import get_versions

__author__ = "Paul Saxe"
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
