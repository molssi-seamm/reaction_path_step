# -*- coding: utf-8 -*-
"""
Control parameters for the Reaction Path step in a SEAMM flowchart
"""

import logging
import seamm

logger = logging.getLogger(__name__)


class ReactionPathParameters(seamm.Parameters):
    """
    The control parameters for Reaction Path.

    You need to replace the "time" entry in dictionary below these comments with the
    definitions of parameters to control this step. The keys are parameters for the
    current plugin,the values are dictionaries as outlined below.

    Examples
    --------
    ::

        parameters = {
            "time": {
                "default": 100.0,
                "kind": "float",
                "default_units": "ps",
                "enumeration": tuple(),
                "format_string": ".1f",
                "description": "Simulation time:",
                "help_text": ("The time to simulate in the dynamics run.")
            },
        }

    parameters : {str: {str: str}}
        A dictionary containing the parameters for the current step.
        Each key of the dictionary is a dictionary that contains the
        the following keys:

    parameters["default"] :
        The default value of the parameter, used to reset it.

    parameters["kind"] : enum()
        Specifies the kind of a variable. One of  "integer", "float", "string",
        "boolean", or "enum"

        While the "kind" of a variable might be a numeric value, it may still have
        enumerated custom values meaningful to the user. For instance, if the parameter
        is a convergence criterion for an optimizer, custom values like "normal",
        "precise", etc, might be adequate. In addition, any parameter can be set to a
        variable of expression, indicated by having "$" as the first character in the
        field. For example, $OPTIMIZER_CONV.

    parameters["default_units"] : str
        The default units, used for resetting the value.

    parameters["enumeration"] : tuple
        A tuple of enumerated values.

    parameters["format_string"] : str
        A format string for "pretty" output.

    parameters["description"] : str
        A short string used as a prompt in the GUI.

    parameters["help_text"] : str
        A longer string to display as help for the user.

    See Also
    --------
    ReactionPath, TkReactionPath, ReactionPath ReactionPathParameters, ReactionPathStep
    """

    parameters = {
        "approach": {
            "default": "Nudged Elastic Band (NEB)",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "Interpolate path",
                "Nudged Elastic Band (NEB)",
            ),
            "format_string": "",
            "description": "Approach:",
            "help_text": "The approach or method for determining the reaction path.",
        },
        "interpolation method": {
            "default": "",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "Image Dependent Pair Potential (IDPP)",
                "Linear",
            ),
            "format_string": "",
            "description": "Interpolation method:",
            "help_text": "How to interpolate needed structures at the beginning",
        },
        "remove rotation and translation": {
            "default": "once before starting",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("once before starting", "every step", "no"),
            "format_string": "",
            "description": "Remove rotation and translation:",
            "help_text": (
                "Whether to remove any rotation and translation between reactants and "
                "products"
            ),
        },
        "number of intermediate structures": {
            "default": 5,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Number of intermediates:",
            "help_text": (
                "The number of intermediate structures between reactants and products."
            ),
        },
        "neb method": {
            "default": "NEB",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("NEB", "AutoNEB"),
            "format_string": "",
            "description": "Method:",
            "help_text": "The NEB method to use.",
        },
        "neb algorithm": {
            "default": "aseneb",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "ase neb",
                "improved tangent",
                "eb",
                "spline",
                "string",
            ),
            "format_string": "",
            "description": "Algorithm:",
            "help_text": "The NEB algorithm to use.",
        },
        "reactant": {
            "default": "",
            "kind": "string",
            "default_units": "",
            "enumeration": ("current configuration", "<system>/<configuration>"),
            "format_string": "",
            "description": "Reactant:",
            "help_text": "The reactant structure.",
        },
        "product": {
            "default": "",
            "kind": "string",
            "default_units": "",
            "enumeration": ("current configuration", "<system>/<configuration>"),
            "format_string": "",
            "description": "Product:",
            "help_text": "The product structure.",
        },
        "intermediate structures": {
            "default": "interpolation",
            "kind": "string",
            "default_units": "",
            "enumeration": (
                "interpolatation",
                "intermediate*",
                "<system>/<configuration>, ...",
            ),
            "format_string": "",
            "description": "Intermediates:",
            "help_text": "The intermediate structures.",
        },
        "spring constant": {
            "default": 500,
            "kind": "float",
            "default_units": "kJ/mol/Å^2",
            "enumeration": tuple(),
            "format_string": ".g",
            "description": "Spring constant:",
            "help_text": "The force constant of the springs between images in the NEB.",
        },
        "number of active images": {
            "default": 3,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Number of active images:",
            "help_text": "The number of simultaneously active images to use.",
        },
        "climbing image": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Use climbing image (CI-NEB):",
            "help_text": "Whether to use a climbing image after convergence.",
        },
        "max climbing steps": {
            "default": "100",
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Maximum # of CI-NEB steps:",
            "help_text": (
                "The maximum number of steps to take in the climbing image phase."
            ),
        },
        "neb optimizer": {
            "default": "BFGS",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "BFGS",
                "FIRE",
                "MDMin",
                "GPMin",
            ),
            "format_string": "",
            "description": "NEB optimizer:",
            "help_text": "The optimizer to use in the NEB part.",
        },
        "initial optimizer": {
            "default": "BFGSLineSearch",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "BFGSLineSearch",
                "BFGS",
                "FIRE",
            ),
            "format_string": "",
            "description": "Reactant/product optimizer:",
            "help_text": "The optimizer to use for the reactant and product.",
        },
        "convergence": {
            "default": 100.0,
            "kind": "float",
            "default_units": "kJ/mol/Å",
            "enumeration": tuple(),
            "format_string": ".g",
            "description": "Convergence criterion:",
            "help_text": "The criterion for convergence of the optimizer.",
        },
        "max steps": {
            "default": "300",
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Maximum # of steps:",
            "help_text": "The maximum number of steps to take.",
        },
        "continue if not converged": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Continue if not converged:",
            "help_text": "Whether to stop if the optimizer does not converge.",
        },
        "on success": {
            "default": "keep last subdirectory",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "keep last subdirectory",
                "keep all subdirectories",
                "delete all subdirectories",
            ),
            "format_string": "",
            "description": "On success:",
            "help_text": "Which subdirectories to keep.",
        },
        "on error": {
            "default": "keep all subdirectories",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "keep last subdirectory",
                "keep all subdirectories",
                "delete all subdirectories",
            ),
            "format_string": "",
            "description": "On error:",
            "help_text": "Which subdirectories to keep if there is an error.",
        },
        "results": {
            "default": {},
            "kind": "dictionary",
            "default_units": None,
            "enumeration": tuple(),
            "format_string": "",
            "description": "results",
            "help_text": "The results to save to variables or in tables.",
        },
    }

    def __init__(self, defaults={}, data=None):
        """
        Initialize the parameters, by default with the parameters defined above

        Parameters
        ----------
        defaults: dict
            A dictionary of parameters to initialize. The parameters
            above are used first and any given will override/add to them.
        data: dict
            A dictionary of keys and a subdictionary with value and units
            for updating the current, default values.

        Returns
        -------
        None
        """

        logger.debug("ReactionPathParameters.__init__")

        super().__init__(
            defaults={
                **ReactionPathParameters.parameters,
                **seamm.standard_parameters.structure_handling_parameters,
                **defaults,
            },
            data=data,
        )

        # Do any local editing of defaults
        tmp = self["structure handling"]
        tmp.default = "Create a new system and configuration"

        tmp = self["subsequent structure handling"]
        tmp.default = "Create a new configuration"

        tmp = self["system name"]
        tmp._data["enumeration"] = (*tmp.enumeration, "RXN")
        tmp.default = "RXN"

        tmp = self["configuration name"]
        tmp._data["enumeration"] = (*tmp.enumeration, "image name")
        tmp.default = "image name"
