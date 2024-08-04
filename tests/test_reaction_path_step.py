#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `reaction_path_step` package."""

import pytest  # noqa: F401
import reaction_path_step  # noqa: F401


def test_construction():
    """Just create an object and test its type."""
    result = reaction_path_step.ReactionPath()
    assert (
        str(type(result)) == "<class 'reaction_path_step.reaction_path.ReactionPath'>"
    )
