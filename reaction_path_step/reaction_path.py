# -*- coding: utf-8 -*-

"""Non-graphical part of the Reaction Path step in a SEAMM flowchart
"""

import logging
from pathlib import Path
import pkg_resources
import pprint  # noqa: F401
import sys

import reaction_path_step
import molsystem
import seamm
from seamm_util import ureg, Q_  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

# In addition to the normal logger, two logger-like printing facilities are
# defined: "job" and "printer". "job" send output to the main job.out file for
# the job, and should be used very sparingly, typically to echo what this step
# will do in the initial summary of the job.
#
# "printer" sends output to the file "step.out" in this steps working
# directory, and is used for all normal output from this step.

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("Reaction Path")

# Add this module's properties to the standard properties
path = Path(pkg_resources.resource_filename(__name__, "data/"))
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)


class ReactionPath(seamm.Node):
    """
    The non-graphical part of a Reaction Path step in a flowchart.

    Attributes
    ----------
    parser : configargparse.ArgParser
        The parser object.

    options : tuple
        It contains a two item tuple containing the populated namespace and the
        list of remaining argument strings.

    subflowchart : seamm.Flowchart
        A SEAMM Flowchart object that represents a subflowchart, if needed.

    parameters : ReactionPathParameters
        The control parameters for Reaction Path.

    See Also
    --------
    TkReactionPath,
    ReactionPath, ReactionPathParameters
    """

    def __init__(
        self,
        flowchart=None,
        title="Reaction Path",
        namespace="org.molssi.seamm.reaction_path",
        extension=None,
        logger=logger,
    ):
        """A step for Reaction Path in a SEAMM flowchart.

        You may wish to change the title above, which is the string displayed
        in the box representing the step in the flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.

        title: str
            The name displayed in the flowchart.
        namespace : str
            The namespace for the plug-ins of the subflowchart
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes

        Returns
        -------
        None
        """
        logger.debug(f"Creating Reaction Path {self}")
        self.subflowchart = seamm.Flowchart(
            parent=self, name="Reaction Path", namespace=namespace
        )  # yapf: disable

        super().__init__(
            flowchart=flowchart,
            title="Reaction Path",
            extension=extension,
            module=__name__,
            logger=logger,
        )  # yapf: disable

        self._metadata = reaction_path_step.metadata

    @property
    def version(self):
        """The semantic version of this module."""
        return reaction_path_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return reaction_path_step.__git_revision__

    def set_id(self, node_id):
        """Set the id for node to a given tuple"""
        self._id = node_id

        # and set our subnodes
        self.subflowchart.set_ids(self._id)

        return self.next()

    def description_text(self, P=None):
        """Create the text description of what this step will do.
        The dictionary of control values is passed in as P so that
        the code can test values, etc.

        Parameters
        ----------
        P: dict
            An optional dictionary of the current values of the control
            parameters.
        Returns
        -------
        str
            A description of the current step.
        """
        self.subflowchart.root_directory = self.flowchart.root_directory

        # Get the first real node
        node = self.subflowchart.get_node("1").next()

        text = self.header + "\n\n"
        while node is not None:
            try:
                text += __(node.description_text(), indent=3 * " ").__str__()
            except Exception as e:
                print(f"Error describing reaction_path flowchart: {e} in {node}")
                logger.critical(
                    f"Error describing reaction_path flowchart: {e} in {node}"
                )
                raise
            except:  # noqa: E722
                print(
                    "Unexpected error describing reaction path flowchart: "
                    f"{sys.exc_info()[0]} in {str(node)}"
                )
                self.logger.critical(
                    "Unexpected error describing reaction path flowchart: "
                    f"{sys.exc_info()[0]} in {str(node)}"
                )
                raise
            text += "\n"
            node = node.next()

        return text

    def run(self):
        """Run a Reaction Path step.

        Parameters
        ----------
        None

        Returns
        -------
        seamm.Node
            The next node object in the flowchart.
        """
        next_node = super().run(printer)
        # Get the first real node
        node = self.subflowchart.get_node("1").next()

        input_data = []
        while node is not None:
            keywords = node.get_input()
            input_data.append(" ".join(keywords))
            node = node.next()

        files = {"molssi.dat": "\n".join(input_data)}
        logger.info("molssi.dat:\n" + files["molssi.dat"])

        local = seamm.ExecLocal()
        result = local.run(
            cmd=["Reaction Path", "-in", "molssi.dat"], files=files, return_files=[]
        )  # yapf: disable

        if result is None:
            logger.error("There was an error running Reaction Path")
            return None

        logger.debug("\n" + pprint.pformat(result))

        logger.info("stdout:\n" + result["stdout"])
        if result["stderr"] != "":
            logger.warning("stderr:\n" + result["stderr"])

        # Analyze the results
        self.analyze()
        # Add other citations here or in the appropriate place in the code.
        # Add the bibtex to data/references.bib, and add a self.reference.cite
        # similar to the above to actually add the citation to the references.

        return next_node

    def analyze(self, indent="", **kwargs):
        """Do any analysis of the output from this step.

        Also print important results to the local step.out file using
        "printer".

        Parameters
        ----------
        indent: str
            An extra indentation for the output
        """
        # Get the first real node
        node = self.subflowchart.get_node("1").next()

        # Loop over the subnodes, asking them to do their analysis
        while node is not None:
            for value in node.description:
                printer.important(value)
                printer.important(" ")

            node.analyze()

            node = node.next()
