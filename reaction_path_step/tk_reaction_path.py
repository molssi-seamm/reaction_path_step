# -*- coding: utf-8 -*-

"""The graphical part of a Reaction Path step"""

import pprint  # noqa: F401
import tkinter as tk  # noqa: F401
import tkinter.ttk as ttk

from .reaction_path_parameters import ReactionPathParameters
import seamm
from seamm_util import ureg, Q_, units_class  # noqa: F401
import seamm_widgets as sw


class TkReactionPath(seamm.TkNode):
    """
    The graphical part of a Reaction Path step in a flowchart.

    Attributes
    ----------
    tk_flowchart : TkFlowchart = None
        The flowchart that we belong to.
    node : Node = None
        The corresponding node of the non-graphical flowchart
    canvas: tkCanvas = None
        The Tk Canvas to draw on
    dialog : Dialog
        The Pmw dialog object
    x : int = None
        The x-coordinate of the center of the picture of the node
    y : int = None
        The y-coordinate of the center of the picture of the node
    w : int = 200
        The width in pixels of the picture of the node
    h : int = 50
        The height in pixels of the picture of the node
    self[widget] : dict
        A dictionary of tk widgets built using the information
        contained in reaction_path_parameters.py

    See Also
    --------
    ReactionPath, TkReactionPath,
    ReactionPathParameters,
    """

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        namespace="org.molssi.seamm.tk",
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50,
    ):
        """
        Initialize a graphical node.

        Parameters
        ----------
        tk_flowchart: Tk_Flowchart
            The graphical flowchart that we are in.
        node: Node
            The non-graphical node for this step.
        namespace: str
            The stevedore namespace for finding sub-nodes.
        canvas: Canvas
           The Tk canvas to draw on.
        x: float
            The x position of the nodes center on the canvas.
        y: float
            The y position of the nodes cetner on the canvas.
        w: float
            The nodes graphical width, in pixels.
        h: float
            The nodes graphical height, in pixels.

        Returns
        -------
        None
        """
        self.namespace = namespace
        self.dialog = None

        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h,
        )
        self.create_dialog()

    def create_dialog(self):
        """
        Create the dialog. A set of widgets will be chosen by default
        based on what is specified in the ReactionPathParameters
        module.

        Parameters
        ----------
        None

        Returns
        -------
        None

        See Also
        --------
        TkReactionPath.reset_dialog
        """

        frame = super().create_dialog(title="Reaction Path", widget="notebook")
        # make it large!
        screen_w = self.dialog.winfo_screenwidth()
        screen_h = self.dialog.winfo_screenheight()
        w = int(0.9 * screen_w)
        h = int(0.8 * screen_h)
        x = int(0.05 * screen_w / 2)
        y = int(0.1 * screen_h / 2)

        self.dialog.geometry(f"{w}x{h}+{x}+{y}")

        # Add a frame for the flowchart
        notebook = self["notebook"]
        flowchart_frame = ttk.Frame(notebook)
        self["flowchart frame"] = flowchart_frame
        notebook.add(flowchart_frame, text="Flowchart", sticky=tk.NSEW)

        self.tk_subflowchart = seamm.TkFlowchart(
            master=flowchart_frame,
            flowchart=self.node.subflowchart,
            namespace=self.namespace,
        )
        self.tk_subflowchart.draw()

        # Fill in the control parameters
        # Shortcut for parameters
        P = self.node.parameters

        # reaction path frame to isolate widgets
        frame = self["reaction path frame"] = ttk.LabelFrame(
            self["frame"],
            borderwidth=4,
            relief="sunken",
            text="Reaction Path Parameters",
            labelanchor="n",
            padding=10,
        )

        for key in ReactionPathParameters.parameters:
            if key not in ("results",):
                self[key] = P[key].widget(frame)

        # and binding to change as needed
        for key in ("approach", "neb method"):
            self[key].combobox.bind(
                "<<ComboboxSelected>>", self.reset_reaction_path_frame
            )

        # structure handling frame to isolate widgets
        frame = self["handling frame"] = ttk.LabelFrame(
            self["frame"],
            borderwidth=4,
            relief="sunken",
            text="Structure Handling",
            labelanchor="n",
            padding=10,
        )

        widgets = []
        row = 0
        for key in (
            "structure handling",
            "subsequent structure handling",
            "system name",
            "configuration name",
        ):
            self[key] = P[key].widget(frame)
            self[key].grid(row=row, column=0, sticky=tk.EW)
            widgets.append(self[key])
            row += 1
        sw.align_labels(widgets, sticky=tk.E)

        # and lay them out
        self.reset_dialog()

        self.setup_results()

    def reset_dialog(self, widget=None):
        """Layout the widgets in the dialog.

        The widgets are chosen by default from the information in
        Diffusivity parameters.

        This function simply lays them out row by row with
        aligned labels. You may wish a more complicated layout that
        is controlled by values of some of the control parameters.
        If so, edit or override this method

        Parameters
        ----------
        widget : Tk Widget = None

        Returns
        -------
        None

        See Also
        --------
        TkDiffusivity.create_dialog
        """

        # Remove any widgets previously packed
        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0

        self["reaction path frame"].grid(row=row, column=0, sticky=tk.EW, pady=10)
        row += 1
        self.reset_reaction_path_frame()

        self["handling frame"].grid(row=row, column=0, sticky=tk.EW, pady=10)
        row += 1

        return row

    def reset_reaction_path_frame(self, widget=None):
        """Layout the widgets in the reaction path frame
        as needed for the current state"""

        approach = self["approach"].get()
        remove = self["remove rotation and translation"].get()

        frame = self["reaction path frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0

        # Main controls
        widgets = []
        widgets2 = []
        for key in ("approach",):
            self[key].grid(row=row, column=0, columnspan=2, sticky=tk.W)
            widgets.append(self[key])
            row += 1

        if approach == "Interpolate path":
            self["remove rotation and translation"].config(values=("yes", "no"))
            if not self.is_expr(remove) and remove != "no":
                self["remove rotation and translation"].set("yes")

            for key in (
                "reactant",
                "product",
                "interpolation method",
                "number of intermediate structures",
                "remove rotation and translation",
            ):
                self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
                widgets2.append(self[key])
                row += 1
        elif approach == "Nudged Elastic Band":
            values = ("once before starting", "every step", "no")
            self["remove rotation and translation"].config(values=values)
            if not self.is_expr(remove) and remove not in values:
                self["remove rotation and translation"].set(values[0])

            method = self["neb method"].get().lower()
            if "auto" in method:
                for key in (
                    "neb method",
                    "reactant",
                    "product",
                    "interpolation method",
                    "number of active images",
                    "number of intermediate structures",
                    "convergence",
                    "max steps",
                    "climbing image",
                    "max climbing steps",
                    "initial optimizer",
                    "continue if not converged",
                    "on success",
                    "on error",
                ):
                    self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
                    widgets2.append(self[key])
                    row += 1
            else:
                for key in (
                    "neb method",
                    "neb algorithm",
                    "reactant",
                    "product",
                    "intermediate structures",
                    "interpolation method",
                    "number of intermediate structures",
                    "remove rotation and translation",
                    "max steps",
                    "convergence",
                    "climbing image",
                    "max climbing steps",
                    "neb optimizer",
                    "spring constant",
                    "continue if not converged",
                    "on success",
                    "on error",
                ):
                    self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
                    widgets2.append(self[key])
                    row += 1

        w1 = sw.align_labels(widgets, sticky=tk.E)
        w2 = sw.align_labels(widgets2, sticky=tk.E)
        frame.columnconfigure(0, minsize=w1 - w2 + 50)

    def right_click(self, event):
        """
        Handles the right click event on the node.

        Parameters
        ----------
        event : Tk Event

        Returns
        -------
        None

        See Also
        --------
        TkReactionPath.edit
        """

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)
