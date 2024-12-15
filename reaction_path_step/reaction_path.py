# -*- coding: utf-8 -*-

"""Non-graphical part of the Reaction Path step in a SEAMM flowchart
"""

import contextlib
from datetime import datetime
import json
import logging
from pathlib import Path
import pkg_resources
import pprint  # noqa: F401
import shutil
import sys
import time
import traceback

from ase import Atoms as ASE_Atoms
from ase.build import (
    minimize_rotation_and_translation as ASE_minimize_rotation_and_translation,
)
from ase.calculators.calculator import all_changes as ASE_all_changes
from ase.mep import (
    AutoNEB as ASE_AutoNEB,
    interpolate as ASE_interpolate,
    idpp_interpolate as ASE_idpp_interpolate,
    NEB as ASE_NEB,
)
import ase.optimize as ASE_Optimize
import numpy as np
from tabulate import tabulate

import reaction_path_step
import molsystem
import read_structure_step
import seamm
from seamm_ase import SEAMM_Calculator
from seamm_util import Q_, units_class, getParser
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


class OutputHandler:
    def __init__(self, *args, keep_newlines=False):
        """Capture stdout and pass to callbacks.

        This is useful for modules or libraries that print information
        directly.
        """
        self.redirector = contextlib.redirect_stdout(self)
        self.callbacks = [*args]
        self.keep_newlines = keep_newlines
        self._buffer = ""

    def add_callback(self, cb):
        """Add a new callback."""
        self.callbacks.append(cb)

    def callbacks(self):
        """Return a list of all the callbacks."""
        return self.callbacks

    def remove_callback(self, cb):
        """Remove a callback."""
        self.callbacks.remove(cb)

    def write(self, txt):
        """Pass the text to the callbacks.

        Buffer the text until it ends in a newline!
        """
        if txt[-1] == "\n":
            # End of line!
            if self.keep_newlines:
                self._buffer += txt
            else:
                self._buffer += txt[0:-1]
            for cb in self.callbacks:
                cb(self._buffer)
            self._buffer = ""
        else:
            self._buffer += txt

    def flush(self):
        pass

    def __enter__(self):
        self.redirector.__enter__()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        # let contextlib do any exception handling here
        self.redirector.__exit__(exc_type, exc_value, traceback)


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
        namespace="org.molssi.seamm",
        extension=None,
        logger=logger,
    ):
        """A step for Reaction Path in a SEAMM flowchart.

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
        logger.setLevel(logging.DEBUG)
        logger.debug(f"Creating Reaction Path {self}")
        self.subflowchart = seamm.Flowchart(
            parent=self, name="Reaction Path", namespace=namespace
        )

        super().__init__(
            flowchart=flowchart,
            title="Reaction Path",
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._metadata = reaction_path_step.metadata
        self.parameters = reaction_path_step.ReactionPathParameters()
        self.optimizer = None
        self._phase = "starting"
        self._images = {}
        self._step = 0
        self._file_handler = None
        self._working_configuration = None
        self._working_data = {}
        self._working_directory = None
        self._data = {}
        self._results = {}
        self._logfile = None
        self._outputfile = None

        self._current_energy = {}
        self._last_coordinates = {}
        self._last_gradients = {}
        self._last_step = None
        self._paths = {}

    @property
    def version(self):
        """The semantic version of this module."""
        return reaction_path_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return reaction_path_step.__git_revision__

    def analyze(self, indent="", step=None, **kwargs):
        """Do any analysis of the output from this step.

        Also print important results to the local step.out file using
        "printer".

        Parameters
        ----------
        indent: str
            An extra indentation for the output
        """
        # Get the energies and path as well as the fit
        structures, path, E, fit_path, fit_E, lines = self.energy_profile(step)

        # Plot to the iteration (step) directory
        self.plot_path(
            self._last_step,
            self._working_directory / f"RXN_path_{self._last_step:04d}.graph",
            path,
            E,
            None,
            None,
            lines,
        )
        # and to the main directory
        self.plot_path(
            self._last_step,
            Path(self.directory) / "RXN_path.graph",
            path,
            E,
            None,
            None,
            lines,
        )

        # Find any maxima along the path
        tol = 1
        Elast = E[0]
        Emin = Elast
        Smin = structures[0]
        Emax = Elast
        Smax = structures[1]
        direction = 0
        extrema = {}
        for structure, energy in zip(structures[1:], E[1:]):
            if energy > Elast + tol:
                if direction <= 0:
                    extrema[Smin] = {"type": "minimum", "energy": Emin}
                    direction = 1
                Emax = energy
                Smax = structure
            elif energy < Elast - tol:
                if direction >= 0:
                    extrema[Smax] = {"type": "maximum", "energy": Emax}
                    direction = -1
                Emin = energy
                Smin = structure
            Elast = energy
        # Last point
        if direction > 0:
            extrema[Smax] = {"type": "maximum", "energy": Emax}
        else:
            extrema[Smin] = {"type": "minimum", "energy": Emin}

        if len(self._paths) == 0:
            self._paths["Step"] = []
            for structure in structures:
                self._paths[structure] = []

        self._paths["Step"].append(step)
        for structure, energy in zip(structures, E):
            if structure in extrema:
                if extrema[structure]["type"] == "maximum":
                    key = "+"
                else:
                    key = "-"
            else:
                key = " "
            self._paths[structure].append(f"{energy:.1f} {key}")

        # Write the paths to an output file
        tabulate.PRESERVE_WHITESPACE = True
        tmp = tabulate(
            self._paths,
            headers="keys",
            tablefmt="rounded_outline",
            disable_numparse=True,
            stralign="right",
            numalign="center",
        )
        tabulate.PRESERVE_WHITESPACE = False
        filepath = Path(self.directory) / "RXN_path.txt"
        filepath.write_text(tmp)

        # Write and SDF file with the key structures
        system = self._working_configuration.system
        configurations = []
        names = []
        configuration = system.get_configuration("reactants")
        names.append(configuration.name)
        configurations.append(configuration)

        pt = 0
        for structure, data in extrema.items():
            if structure == "reactants" or structure == "products":
                continue
            configuration = system.get_configuration(structure)
            names.append(configuration.name)
            pt += 1
            if data["type"] == "maximum":
                title = f"TS_{pt}"
            elif data["type"] == "minimum":
                title = f"IM_{pt}"
            configuration.name = title
            configurations.append(configuration)

        configuration = system.get_configuration("products")
        names.append(configuration.name)
        configurations.append(configuration)

        filepath = Path(self.directory) / "RXN_path.sdf"
        self._write_sdf(filepath, configurations)

        filepath = Path(self.directory) / "neb" / f"RXN_path_{step:04d}.sdf"
        self._write_sdf(filepath, configurations)

        # Put the configuration names back
        for configuration, name in zip(configurations, names):
            configuration.name = name

    def attach_calculators(self, images):
        """Attach calculators to the images in the reaction path."""
        printer.job(f"phase {self._phase}, {images=}")
        if "iteration" in self._phase:
            printer.job(f"images are {self._working_data['images']}")
            for i, image in enumerate(images):
                image.calc = SEAMM_Calculator(self)
                image_no = self._working_data["images"][i + 1]
                if image_no in self._images:
                    if self._images[image_no] != image.calc:
                        raise RuntimeError(f"Error in image no {image_no}")
                else:
                    self._images[image_no] = image.calc
                printer.job(f"   image {image_no} --> {self._images[image_no]}")
        else:
            for i in range(len(images)):
                images[i].calc = SEAMM_Calculator(self)

    def auto_neb(self, P):
        """Do the NEB calculation.

        Parameters
        ----------
        P : dict()
            The control parameters for the step.
        """
        self._data = {
            "step": [],
            "energy": [],
            "max_force": [],
            "rms_force": [],
            "max_step": [],
        }
        self._working_data = {}
        self._last_coordinates = None
        self._step = 0

        wd = Path(self.directory)

        # Optimize the reaction path
        if P["neb method"] == "AutoNEB":
            # Get the reactants and products configurations
            reactants, products = self.get_reactants_and_products(P)

            # The default maximum number of steps may depend on the number of atoms
            n_atoms = reactants.n_atoms
            max_steps = P["max steps"]
            if isinstance(max_steps, str) and "natoms" in max_steps:
                tmp = max_steps.split()
                if "natoms" in tmp[0]:
                    max_steps = int(tmp[1]) * n_atoms
                else:
                    max_steps = int(tmp[0]) * n_atoms

            # Write the initial structures to disk

            indent = self.indent + 4 * " "
            printer.normal(__("    Optimizing reactants structure", indent=indent))

            self._phase = "reactants optimization"
            self._logfile = wd / "reactants_optimization.log"
            self._outputfile = wd / "reactants_optimization.out"
            self._working_configuration = reactants
            self._working_directory = wd / "reactants"
            self._step = 0

            # Set the default system and configuration to the reactants
            system_db = reactants.system_db

            system_db.system = reactants.system
            system_db.system.configuration = reactants

            ASE_reactants = ASE_Atoms(
                "".join(reactants.atoms.symbols), positions=reactants.atoms.coordinates
            )
            ASE_reactants.calc = SEAMM_Calculator(self)
            with OutputHandler(self.log_calculator):
                with ASE_Optimize.BFGS(
                    ASE_reactants, trajectory=str(wd / "neb000.traj")
                ) as opt:
                    opt.run(fmax=P["convergence"].m_as("eV/Å"), steps=max_steps)
            reactants.atoms.set_coordinates(ASE_reactants.positions)
            printer.normal(__("       updated the reactants structure", indent=indent))
            self.write_structure(reactants, wd / "reactants.mmcif")

            printer.normal(__("Optimizing products structure", indent=indent))

            self._phase = "products optimization"
            self._logfile = wd / "products_optimization.log"
            self._outputfile = wd / "products_optimization.out"
            self._working_configuration = products
            self._working_directory = wd / "products"
            self._step = 0

            # Set the default system and configuration to the products
            system_db.system = products.system
            system_db.system.configuration = products

            ASE_products = ASE_Atoms(
                "".join(products.atoms.symbols), positions=products.atoms.coordinates
            )
            ASE_products.calc = SEAMM_Calculator(self)
            with OutputHandler(self.log_calculator):
                with ASE_Optimize.BFGS(
                    ASE_products, trajectory=str(wd / "neb001.traj")
                ) as opt:
                    opt.run(fmax=P["convergence"].m_as("eV/Å"), steps=max_steps)
            products.atoms.set_coordinates(ASE_products.positions)
            printer.normal(__("       updated the products structure", indent=indent))
            self.write_structure(products, wd / "products.mmcif")

            # Make a new configuration to handle the images.
            _, configuration = self.get_system_configuration(P)
            self._working_configuration = configuration

            printer.normal(__("Starting the AutoNEB calculation", indent=indent))

            self._phase = "AutoNEB"
            self._logfile = wd / "AutoNEB.log"
            self._outputfile = wd / "AutoNEB.out"
            self._working_directory = wd / "neb"
            self._step = 0

            # At the moment only two choices supported for the interpolation
            if "IDPP" in P["interpolation method"].lower():
                interpolation_method = "idpp"
            else:
                interpolation_method = "linear"

            driver = ASE_AutoNEB(
                self.attach_calculators,
                prefix=str(wd / "neb"),
                method="eb",  # "aseneb", "improvedtangent", "eb"
                interpolate_method=interpolation_method,
                n_simul=P["number of active images"],
                n_max=P["number of intermediate structures"] + 2,
                climb=P["climbing image"],
                fmax=P["convergence"].m_as("eV/Å"),
                maxsteps=[max_steps, P["max climbing steps"]],
                k=0.5,  # spring constant along NEB path
                optimizer=ASE_Optimize.BFGS,
                parallel=False,
            )

            # Set up the NEB calculation
            raise_exception = False
            error = None
            tic = time.perf_counter_ns()
            try:
                with OutputHandler(self.log_calculator, self.auto_neb_handler):
                    ASE_images = driver.run()
            except Exception as err:  # noqa: F841
                self.logger.error("Caught exception: ", err)
                error = err
                raise_exception = True
            else:
                best = ASE_images[0]
                best_e = ASE_images[0].get_total_energy()
                for i, image in enumerate(ASE_images):
                    if image.get_total_energy() > best_e:
                        best = image
                        best_e = image.get_total_energy()
                configuration.atoms.set_coordinates(best.positions)
            finally:
                toc = time.perf_counter_ns()
                self._results["t_elapsed"] = round((toc - tic) * 1.0e-9, 3)

                # Print the results
                # self.analyze()

                # Store results to db, variables, tables, and json as requested
                self.store_results(
                    configuration=self._working_configuration,
                    data=self._results,
                )

            # Clean up the subdirectories
            if raise_exception:
                keep = P["on error"]
                if keep == "delete all subdirectories":
                    subdirectories = wd.glob("step_*")
                    for subdirectory in subdirectories:
                        shutil.rmtree(subdirectory)
                elif keep == "keep last subdirectory":
                    subdirectories = wd.glob("step_*")
                    subdirectories = sorted(subdirectories)
                    for subdirectory in subdirectories[:-1]:
                        shutil.rmtree(subdirectory)
                raise error from None
            else:
                keep = P["on success"]
                if keep == "delete all subdirectories":
                    subdirectories = wd.glob("step_*")
                    for subdirectory in subdirectories:
                        shutil.rmtree(subdirectory)
                elif keep == "keep last subdirectory":
                    subdirectories = wd.glob("step_*")
                    subdirectories = sorted(subdirectories)
                    for subdirectory in subdirectories[:-1]:
                        shutil.rmtree(subdirectory)
        else:
            raise ValueError(f"Unknown NEB method: {P['neb method']}")

    def auto_neb_handler(self, txt):
        """Track what is going on in AutoNEB from the output"""
        if txt.startswith("Start of evaluation of the initial images"):
            self._phase = "initial images"

        if txt.startswith("Now starting iteration"):
            it = int(txt.split()[3])
            tmp = txt.split("[")[1][0:-1].split(", ")
            n_images = len(tmp)
            self._working_data["images"] = [int(i) for i in tmp]
            self._phase = f"iteration {it}"
            self._working_directory = Path(self.directory) / "neb" / f"Iter_{it:02d}"
            printer.normal(
                __(
                    f"    Iteration {it} using {n_images} images",
                    indent=self.indent + 4 * " ",
                )
            )

    def calculator_for_neb(
        self,
        calculator,
        properties=["energy"],
        system_changes=ASE_all_changes,
    ):
        """Create a calculator for the ASE NEB method

        Parameters
        ----------
        ase : ase.calculators.calculator.Calculator
            The ASE calculator we are working for
        properties : list of str
            The properties to calculate.
        system_changes : int
            The changes to the system.

        Returns
        -------
        results : dict
            The dictionary of results from the calculation.
        """
        fmt = "04d"

        wd = self._working_directory
        wd.mkdir(parents=True, exist_ok=True)

        # Get the configuration that we are working on
        self._working_configuration = configuration = calculator.configuration

        # And set as current
        configuration.system_db.system = system = configuration.system
        system.configuration = configuration

        # Get the name to use for files, etc.
        name = calculator.name.replace(" ", "_")
        self._data["structure"].append(name)

        if self.optimizer is None:
            self._step = 0
        else:
            self._step = self.optimizer.get_number_of_steps()
        self._results["nsteps"] = self._step

        # Create the directory for this step and structure
        step_id = f"step_{self._step:{fmt}}"
        step_dir = self._working_directory / step_id
        calc_dir = step_dir / name
        calc_dir.mkdir(parents=True, exist_ok=True)

        if self._step == self._last_step:
            self._data["step"].append("")
        else:
            self._data["step"].append(self._step)
            if self._last_step is not None:
                self.analyze(step=self._last_step)
            self._last_step = self._step

        # Write so that we can see it
        # self.write_structure(configuration, wd.parent / name)

        calculator.results = {}

        n_atoms = len(calculator.atoms)
        self.logger.debug(f"{n_atoms} atoms in the reaction path")
        positions = calculator.atoms.positions
        self.logger.debug(f"Positions: {positions}")
        cell = calculator.atoms.cell
        self.logger.debug(f"Cell: {cell}")

        # Set the coordinates in the configuration
        configuration.atoms.set_coordinates(positions, fractionals=False)

        # Find the handler for job.out and set the level up
        job_handler = None
        out_handler = None
        for handler in job.handlers:
            if (
                isinstance(handler, logging.FileHandler)
                and "job.out" in handler.baseFilename
            ):
                job_handler = handler
                job_level = job_handler.level
                job_handler.setLevel(printing.JOB)
            elif isinstance(handler, logging.StreamHandler):
                out_handler = handler
                out_level = out_handler.level
                out_handler.setLevel(printing.JOB)

        # Get the first real node
        first_node = self.subflowchart.get_node("1").next()

        # Ensure the nodes have their options
        node = first_node
        while node is not None:
            node.all_options = self.all_options
            node = node.next()

        # And the subflowchart has the executor
        self.subflowchart.executor = self.flowchart.executor

        # A handler for the file
        if self._file_handler is not None:
            self._file_handler.close()
            job.removeHandler(self._file_handler)
        path = calc_dir / "Step.out"
        path.unlink(missing_ok=True)
        self._file_handler = logging.FileHandler(path)
        self._file_handler.setLevel(printing.NORMAL)
        formatter = logging.Formatter(fmt="{message:s}", style="{")
        self._file_handler.setFormatter(formatter)
        job.addHandler(self._file_handler)

        # Add the step to the ids so the directory structure is reasonable
        self.subflowchart.reset_visited()
        tmp = self._working_directory.relative_to(self.directory).parts
        self.set_subids((*self._id, *tmp, step_id, name))

        # Run through the steps in the loop body
        node = first_node
        try:
            while node is not None:
                node = node.run()
        except DeprecationWarning as e:
            printer.normal("\nDeprecation warning: " + str(e))
            traceback.print_exc(file=sys.stderr)
            traceback.print_exc(file=sys.stdout)
        except Exception as e:
            printer.job(f"Caught exception in step {self._step}/{name}: {str(e)}")
            with open(calc_dir / "stderr.out", "a") as fd:
                traceback.print_exc(file=fd)
            raise
        self.logger.debug(f"End of step {self._step}/{name}")

        # Write so that we can see it
        # self.write_structure(configuration, wd.parent / name)

        # Remove any redirection of printing.
        if self._file_handler is not None:
            self._file_handler.close()
            job.removeHandler(self._file_handler)
            self._file_handler = None
        if job_handler is not None:
            job_handler.setLevel(job_level)
        if out_handler is not None:
            out_handler.setLevel(out_level)

        # Get the energy and derivatives
        paths = sorted(calc_dir.glob("**/Results.json"))

        if len(paths) == 0:
            raise RuntimeError(
                "There are no energy and gradients in properties.json for step "
                f"{self._step}/{name} in {calc_dir}."
            )
        else:
            # Find the most recent and assume that is the one wanted
            newest_time = None
            for path in paths:
                with path.open() as fd:
                    data = json.load(fd)
                time = datetime.fromisoformat(data["iso time"])
                if newest_time is None:
                    newest = path
                    newest_time = time
                elif time > newest_time:
                    newest_time = time
                    newest = path
            with newest.open() as fd:
                data = json.load(fd)

        if "energy,units" in data:
            units = data["energy,units"]
        else:
            units = "kJ/mol"
        energy = Q_(data["energy"], units)
        self._data["energy"].append(energy.m_as("kJ/mol"))

        self._results["energy"] = energy.m_as("kJ/mol")
        self._current_energy[calculator.name] = energy.m_as("kJ/mol")

        gradients = data["gradients"]

        if self.logger.isEnabledFor(logging.DEBUG):
            self.logger.debug("\ngradients")
            for i in range(n_atoms):
                self.logger.debug(
                    f"   {gradients[i][0]:8.3f} {gradients[i][1]:8.3f} "
                    f"{gradients[i][2]:8.3f}"
                )

        if "gradients,units" in data:
            funits = data["gradients,units"]
        else:
            funits = "kJ/mol/Å"

        # Get the measures of convergence
        max_force = np.max(np.linalg.norm(gradients, axis=1))
        self._data["max_force"].append(max_force)
        self._results["maximum_gradient"] = Q_(max_force, funits).m_as("kJ/mol/Å")
        rms_force = np.sqrt(np.mean(np.linalg.norm(gradients, axis=1) ** 2))
        self._data["rms_force"].append(rms_force)
        self._results["rms_gradient"] = Q_(rms_force, funits).m_as("kJ/mol/Å")

        if self._step > 1:
            step = positions - self._last_coordinates[name]
            max_step = np.max(np.linalg.norm(step, axis=1))
        else:
            max_step = 0.0
        self._data["max_step"].append(max_step)
        self._results["maximum_step"] = max_step
        self._last_coordinates[name] = np.array(positions)

        # Units!
        gradients = np.array(gradients) * Q_(1.0, funits).to("eV/Å").magnitude
        self._last_gradients[name] = (
            np.array(gradients) * Q_(1.0, funits).to("kJ/mol/Å").magnitude
        )

        calculator.results["energy"] = energy.m_as("eV")
        calculator.results["forces"] = -gradients

        # Log the results
        if self._logfile is not None:
            headers = [
                "Step",
                "Structure",
                f"E ({units})",
                f"Fmax ({funits})",
                f"Frms ({funits})",
                "max step (Å)",
            ]
            tmp = tabulate(
                self._data,
                headers=headers,
                tablefmt="rounded_outline",
                disable_numparse=False,
                floatfmt=".3f",
            )
            with open(self._logfile, "w") as fd:
                fd.write(tmp)
                fd.write("\n")

    def calculator(
        self,
        calculator,
        properties=["energy"],
        system_changes=ASE_all_changes,
    ):
        """Create a calculator for the reaction path step.

        Parameters
        ----------
        ase : ase.calculators.calculator.Calculator
            The ASE calculator we are working for
        properties : list of str
            The properties to calculate.
        system_changes : int
            The changes to the system.

        Returns
        -------
        results : dict
            The dictionary of results from the calculation.
        """
        wd = self._working_directory
        wd.mkdir(parents=True, exist_ok=True)

        # Initialize the data when starting or restarting an optimization
        if self._step == 0:
            self._data = {
                "step": [],
                "energy": [],
                "max_force": [],
                "rms_force": [],
                "max_step": [],
            }

        self._step += 1
        self._results["nsteps"] = self._step
        self._data["step"].append(self._step)
        fmt = "05d"

        calculator.results = {}

        n_atoms = len(calculator.atoms)
        self.logger.debug(f"{n_atoms} atoms in the reaction path")
        positions = calculator.atoms.positions
        self.logger.debug(f"Positions: {positions}")
        cell = calculator.atoms.cell
        self.logger.debug(f"Cell: {cell}")

        # Set the coordinates in the configuration
        self._working_configuration.atoms.set_coordinates(positions, fractionals=False)

        # Find the handler for job.out and set the level up
        job_handler = None
        out_handler = None
        for handler in job.handlers:
            if (
                isinstance(handler, logging.FileHandler)
                and "job.out" in handler.baseFilename
            ):
                job_handler = handler
                job_level = job_handler.level
                job_handler.setLevel(printing.JOB)
            elif isinstance(handler, logging.StreamHandler):
                out_handler = handler
                out_level = out_handler.level
                out_handler.setLevel(printing.JOB)

        # Get the first real node
        first_node = self.subflowchart.get_node("1").next()

        # Ensure the nodes have their options
        node = first_node
        while node is not None:
            node.all_options = self.all_options
            node = node.next()

        # And the subflowchart has the executor
        self.subflowchart.executor = self.flowchart.executor

        # Direct most output to iteration.out
        step_id = f"step_{self._step:{fmt}}"
        step_dir = self._working_directory / step_id
        step_dir.mkdir(parents=True, exist_ok=True)

        # A handler for the file
        if self._file_handler is not None:
            self._file_handler.close()
            job.removeHandler(self._file_handler)
        path = step_dir / "Step.out"
        path.unlink(missing_ok=True)
        self._file_handler = logging.FileHandler(path)
        self._file_handler.setLevel(printing.NORMAL)
        formatter = logging.Formatter(fmt="{message:s}", style="{")
        self._file_handler.setFormatter(formatter)
        job.addHandler(self._file_handler)

        # Add the step to the ids so the directory structure is reasonable
        self.subflowchart.reset_visited()
        tmp = self._working_directory.relative_to(self.directory).parts
        self.set_subids((*self._id, *tmp, step_id))

        # Run through the steps in the loop body
        node = first_node
        try:
            while node is not None:
                node = node.run()
        except DeprecationWarning as e:
            printer.normal("\nDeprecation warning: " + str(e))
            traceback.print_exc(file=sys.stderr)
            traceback.print_exc(file=sys.stdout)
        except Exception as e:
            printer.job(f"Caught exception in step {self._step}: {str(e)}")
            with open(step_dir / "stderr.out", "a") as fd:
                traceback.print_exc(file=fd)
            raise
        self.logger.debug(f"End of step {self._step}")

        # Remove any redirection of printing.
        if self._file_handler is not None:
            self._file_handler.close()
            job.removeHandler(self._file_handler)
            self._file_handler = None
        if job_handler is not None:
            job_handler.setLevel(job_level)
        if out_handler is not None:
            out_handler.setLevel(out_level)

        # Get the energy and derivatives
        paths = sorted(step_dir.glob("**/Results.json"))

        if len(paths) == 0:
            raise RuntimeError(
                "There are no energy and gradients in properties.json for step "
                f"{self._step} in {step_dir}."
            )
        else:
            # Find the most recent and assume that is the one wanted
            newest_time = None
            for path in paths:
                with path.open() as fd:
                    data = json.load(fd)
                time = datetime.fromisoformat(data["iso time"])
                if newest_time is None:
                    newest = path
                    newest_time = time
                elif time > newest_time:
                    newest_time = time
                    newest = path
            with newest.open() as fd:
                data = json.load(fd)

        energy = data["energy"]
        if "energy,units" in data:
            units = data["energy,units"]
        else:
            units = "kJ/mol"
        self._data["energy"].append(energy)

        energy *= Q_(1.0, units).to("eV").magnitude
        self._results["energy"] = Q_(energy, "eV").m_as("kJ/mol")

        gradients = data["gradients"]

        if self.logger.isEnabledFor(logging.DEBUG):
            self.logger.debug("\ngradients")
            for i in range(n_atoms):
                self.logger.debug(
                    f"   {gradients[i][0]:8.3f} {gradients[i][1]:8.3f} "
                    f"{gradients[i][2]:8.3f}"
                )

        if "gradients,units" in data:
            funits = data["gradients,units"]
        else:
            funits = "kJ/mol/Å"

        # Get the measures of convergence
        max_force = np.max(np.linalg.norm(gradients, axis=1))
        self._data["max_force"].append(max_force)
        self._results["maximum_gradient"] = Q_(max_force, funits).m_as("kJ/mol/Å")
        rms_force = np.sqrt(np.mean(np.linalg.norm(gradients, axis=1) ** 2))
        self._data["rms_force"].append(rms_force)
        self._results["rms_gradient"] = Q_(rms_force, funits).m_as("kJ/mol/Å")

        if self._step > 1:
            step = positions - self._last_coordinates
            max_step = np.max(np.linalg.norm(step, axis=1))
        else:
            max_step = 0.0
        self._data["max_step"].append(max_step)
        self._results["maximum_step"] = max_step
        self._last_coordinates = np.array(positions)

        # Units!
        gradients = np.array(gradients) * Q_(1.0, funits).to("eV/Å").magnitude

        calculator.results["energy"] = energy
        calculator.results["forces"] = -gradients

        # Log the results
        if self._logfile is not None:
            headers = [
                "Step",
                f"E ({units})",
                f"Fmax ({funits})",
                f"Frms ({funits})",
                "max step (Å)",
            ]
            tmp = tabulate(
                self._data,
                headers=headers,
                tablefmt="rounded_outline",
                disable_numparse=False,
                floatfmt=".3f",
            )
            with open(self._logfile, "w") as fd:
                fd.write(tmp)
                fd.write("\n")

    def create_parser(self):
        """Setup the command-line / config file parser"""
        parser_name = "reaction-path-step"
        parser = getParser()

        # Remember if the parser exists ... this type of step may have been
        # found before
        parser_exists = parser.exists(parser_name)

        # Create the standard options, e.g. log-level
        super().create_parser(name=parser_name)

        if not parser_exists:
            # Any options for diffusivity itself
            parser.add_argument(
                parser_name,
                "--html",
                action="store_true",
                help="whether to write out html files for graphs, etc.",
            )

        # Now need to walk through the steps in the subflowchart...
        self.subflowchart.reset_visited()
        node = self.subflowchart.get_node("1").next()
        while node is not None:
            node = node.create_parser()

        return self.next()

    def description_text(self, P=None, short=False, natoms=None):
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
        if P is None:
            P = self.parameters.values_to_dict()

        approach = P["approach"].lower()
        if "interpolat" in approach:
            text = "The path between reactants '{reactant}' and products '{product}' "
            text += "will be interpolated using the {interpolation method} method "
            text += "to create {number of intermediate structures} intermediate "
            text += "structures. "
            remove = P["remove rotation and translation"]
            if remove == "no":
                text += "The products will not be translated and rotated to best "
                text += "overlay the reactants."
            elif self.is_expr(remove):
                text += "Whether to rotate and translate the products to best overlay "
                text += "the reactants will be determined by "
                text += "{remove rotation and translation}."
            else:
                text += "The products will be rotated and translated to best overlay "
                text += f"the reactants {remove}."
        elif "neb" in approach or "nudge" in approach:
            text = "The reaction path will be explored using nudged elastic band "
            text += "method ({neb method}) using the {neb algorithm}, converging to "
            text += "{convergence} "

            max_steps = P["max steps"]
            if (
                natoms is not None
                and isinstance(max_steps, int)
                and "natoms" in max_steps
            ):
                tmp = max_steps.split()
                if "natoms" in tmp[0]:
                    max_steps = int(tmp[1]) * natoms
                else:
                    max_steps = int(tmp[0]) * natoms
            text += f"with a maximum of {max_steps} steps."

            remove = P["remove rotation and translation"]
            if self.is_expr(remove):
                text += " Whether to rotate and translate the products to best overlay "
                text += "the reactants will be determined by "
                text += "{remove rotation and translation}."
            elif "once" in remove or "every" in remove:
                text += " The products will be rotated and translated to best overlay "
                text += f"the reactants {remove}."
            elif remove == "no":
                text += " The products will not be translated and rotated to best "
                text += "overlay the reactants."
            else:
                raise ValueError(
                    "Don't understand option to remove rotation and translation: "
                    f"'{remove}'"
                )
            stop = P["continue if not converged"]
            if isinstance(stop, bool) and not stop or stop == "no":
                text += " The workflow will continue if the NEB "
                text += "does not converge."

            # Make sure the subflowchart has the data from the parent flowchart
            self.subflowchart.root_directory = self.flowchart.root_directory
            self.subflowchart.executor = self.flowchart.executor
            self.subflowchart.in_jobserver = self.subflowchart.in_jobserver

            if not short:
                # Get the first real node
                node = self.subflowchart.get_node("1").next()
                text += "\n\n"

                # Now walk through the steps in the subflowchart...
                while node is not None:
                    try:
                        text += __(node.description_text(), indent=3 * " ").__str__()
                    except Exception as e:
                        print(
                            f"Error describing reaction path flowchart: {e} in {node}"
                        )
                        self.logger.critical(
                            f"Error describing reaction path flowchart: {e} in {node}"
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

        return self.header + "\n" + __(text, **P, indent=4 * " ").__str__()

    def energy_profile(self, step):
        """Prepare the energy profile for this step"""
        energies = {}
        for _step, structure, energy in zip(
            self._data["step"], self._data["structure"], self._data["energy"]
        ):
            if _step != "" and _step > step:
                break
            energies[structure] = energy

        npts = len(energies)
        fit_energies = np.empty((npts - 1) * 20 + 1)
        fit_path = np.empty((npts - 1) * 20 + 1)

        # Normalize to the energy of the reactants
        if "reactants" in energies:
            E0 = energies["reactants"]
            for structure in energies.keys():
                energies[structure] -= E0

        # Order structures from reactants to products
        n_intermediates = len(energies) - 2
        structures = ["reactants"]
        for i in range(1, n_intermediates + 1):
            key = f"img_{i}"
            if key not in energies:
                raise RuntimeError(f"Missing intermediate structure {i}.")
            structures.append(key)
        structures.append("products")

        E = [energies[structure] for structure in structures]

        path = [0]
        for i in range(npts - 1):
            s0 = structures[i]
            s1 = structures[i + 1]
            dR = self._last_coordinates[s1] - self._last_coordinates[s0]
            path.append(path[i] + np.sqrt((dR**2).sum()))

        # Prepare a fit, following how ASE does this
        lines = []  # tangent lines
        last_slope = None
        for i in range(npts):
            if i == 0:
                direction = (
                    self._last_coordinates[structures[i + 1]]
                    - self._last_coordinates[structures[i]]
                )
                dpath = 0.5 * path[1]
            elif i == npts - 1:
                direction = (
                    self._last_coordinates[structures[-1]]
                    - self._last_coordinates[structures[-2]]
                )
                dpath = 0.5 * (path[-1] - path[-2])
            else:
                direction = (
                    self._last_coordinates[structures[i + 1]]
                    - self._last_coordinates[structures[i - 1]]
                )
                dpath = 0.25 * (path[i + 1] - path[i - 1])

            direction /= np.linalg.norm(direction)
            slope = (self._last_gradients[structures[i]] * direction).sum()
            x = np.linspace(path[i] - dpath, path[i] + dpath, 3)
            y = E[i] + slope * (x - path[i])
            lines.append((x, y))

            if i > 0:
                s0 = path[i - 1]
                s1 = path[i]
                x = np.linspace(s0, s1, 20, endpoint=False)
                c = np.linalg.solve(
                    np.array(
                        [
                            (1, s0, s0**2, s0**3),
                            (1, s1, s1**2, s1**3),
                            (0, 1, 2 * s0, 3 * s0**2),
                            (0, 1, 2 * s1, 3 * s1**2),
                        ]
                    ),
                    np.array([E[i - 1], E[i], last_slope, slope]),
                )
                y = c[0] + x * (c[1] + x * (c[2] + x * c[3]))
                fit_path[(i - 1) * 20 : i * 20] = x
                fit_energies[(i - 1) * 20 : i * 20] = y

            last_slope = slope

        fit_path[-1] = path[-1]
        fit_energies[-1] = E[-1]

        return structures, path, E, fit_path, fit_energies, lines

    def get_reactants_and_products(self, P):
        """Get the reactants and products systems.

        Parameters
        ----------
        P : dict()
            The control parameters for the step
        """
        system_db = self.get_variable("_system_db")
        if P["reactant"] == "current configuration":
            _, reactants = self.get_system_configuration()
        else:
            if "/" in P["reactant"]:
                sysname, confname = P["reactant"].split("/", 1)
                tmp_system = system_db.get_system(sysname)
            else:
                confname = P["reactant"]
                tmp_system, _ = self.get_system_configuration()
            reactants = tmp_system.get_configuration(confname)
        reactants_name = f"{reactants.system.name}.{reactants.name}"

        if P["product"] == "current configuration":
            _, products = self.get_system_configuration()
        else:
            if "/" in P["product"]:
                sysname, confname = P["product"].split("/", 1)
                tmp_system = system_db.get_system(sysname)
            else:
                confname = P["product"]
                tmp_system, _ = self.get_system_configuration()
            products = tmp_system.get_configuration(confname)
        products_name = f"{products.system.name}.{products.name}"

        if reactants.n_atoms != products.n_atoms:
            msg = "The reactants and products must have the same number of atoms!"
            msg += "\n"
            msg += f"  reactants has {reactants.n_atoms} ({reactants_name})"
            msg += "\n"
            msg += f"   products has {products.n_atoms} ({products_name})"
            printer.job(msg)
            printer.job("")
            raise RuntimeError(msg)

        return reactants, products

    def interpolate(self, P):
        """Interpolate the structures for the reaction

        Parameters
        ----------
        P : dict()
            The control parameters for the step.

        Returns
        -------
        configurations : [molsystem._Configuration]
            The list of configurations along the path from reactants to products
        """
        # Get the reactants and products configurations
        reactants, products = self.get_reactants_and_products(P)

        method = P["interpolation method"].lower()
        if method == "linear" or "idpp" in method:
            # Simplest approach: linear interpolation
            ASE_reactants = ASE_Atoms(
                "".join(reactants.atoms.symbols), positions=reactants.atoms.coordinates
            )
            ASE_reactants.calc = SEAMM_Calculator(
                self, name="reactants", configuration=reactants
            )
            ASE_products = ASE_Atoms(
                "".join(products.atoms.symbols), positions=products.atoms.coordinates
            )
            ASE_products.calc = SEAMM_Calculator(
                self, name="products", configuration=products
            )

            n_intermediates = P["number of intermediate structures"]
            images = [ASE_reactants]
            for _ in range(n_intermediates):
                images += [ASE_reactants.copy()]
            images += [ASE_products]

            if P["remove rotation and translation"] != "no":
                ASE_minimize_rotation_and_translation(images[0], images[-1])

            if "idpp" in method:
                # Use the Image Dependent Pair Potential (IDPP) approach (like LST!)
                ASE_idpp_interpolate(images=images, traj=None, log=None)
            else:
                ASE_interpolate(images)

        # Create a new system with the images
        configurations = []
        system, configuration = self.get_system_configuration(P, same_as=reactants)
        configurations.append(configuration)
        if P["configuration name"] == "image name":
            configuration.name = "reactants"
        configuration.atoms.set_coordinates(images[0].positions)

        for i in range(n_intermediates):
            system, configuration = self.get_system_configuration(P, first=False)
            configurations.append(configuration)
            if P["configuration name"] == "image name":
                configuration.name = f"img_{i+1}"
            configuration.atoms.set_coordinates(images[i + 1].positions)

        system, configuration = self.get_system_configuration(
            P, first=False, same_as=products
        )
        configurations.append(configuration)
        if P["configuration name"] == "image name":
            configuration.name = "products"
        configuration.atoms.set_coordinates(images[-1].positions)

        return configurations

    def log_calculator(self, txt):
        """Direct captured output to the correct file.

        Parameters
        ----------
        txt : str
            The text to write.
        """
        with self._outputfile.open(mode="a") as fd:
            fd.write(txt + "\n")

    def neb(self, P):
        """Do the NEB calculation.

        Parameters
        ----------
        P : dict()
            The control parameters for the step.
        """
        # Initialization
        wd = Path(self.directory)
        indent = self.indent + 4 * " "

        self._data = {
            "step": [],
            "structure": [],
            "energy": [],
            "max_force": [],
            "rms_force": [],
            "max_step": [],
        }
        self._paths = {}
        self._current_energy = {}
        self._last_coordinates = {}
        self._last_gradients = {}
        self._last_step = None
        self._logfile = wd / "NEB.log"
        self._outputfile = wd / "NEB.out"
        self._phase = "NEB"
        self._step = 0
        self._working_data = {}
        self._working_directory = wd / "neb"

        # Get the reactants and products configurations
        reactants, products = self.get_reactants_and_products(P)

        # The default maximum number of steps may depend on the number of atoms
        n_atoms = reactants.n_atoms
        max_steps = P["max steps"]
        if isinstance(max_steps, str) and "natoms" in max_steps:
            tmp = max_steps.split()
            if "natoms" in tmp[0]:
                max_steps = int(tmp[1]) * n_atoms
            else:
                max_steps = int(tmp[0]) * n_atoms

        printer.normal(__("Starting the NEB calculation", indent=indent))

        # And the intermediate structures...
        tmp = P["intermediate structures"]
        if "interpolat" in tmp.lower():
            configurations = self.interpolate(P)
        else:
            raise NotImplementedError("Explicit path not implemented yet!")

        # Set up the images and calculators for the ASE NEB code
        images = []
        for i, configuration in enumerate(configurations):
            if i == 0:
                name = "reactants"
            elif i == len(configurations) - 1:
                name = "products"
            else:
                name = f"img_{i}"
            image = ASE_Atoms(
                "".join(configuration.atoms.symbols),
                positions=configuration.atoms.coordinates,
            )
            image.calc = SEAMM_Calculator(
                self,
                calculator=self.calculator_for_neb,
                name=name,
                configuration=configuration,
            )
            images.append(image)

        # Calculate the energy and forces on the reactants and products
        images[0].calc.get_potential_energy(atoms=images[0])
        images[-1].calc.get_potential_energy(atoms=images[-1])

        # Set up the NEB calculation
        remove = "every" in P["remove rotation and translation"]
        spring_constant = P["spring constant"].m_as("eV/Å^2")
        convergence = P["convergence"].m_as("eV/Å")

        raise_exception = False
        error = None
        tic = time.perf_counter_ns()
        converged = False
        try:
            with OutputHandler(self.log_calculator):
                neb = ASE_NEB(
                    images,
                    k=spring_constant,
                    climb=P["climbing image"],
                    remove_rotation_and_translation=remove,
                    method=P["neb algorithm"].replace(" ", "").lower(),
                )
                match P["neb optimizer"].lower():
                    case "bfgs":
                        optimizer = ASE_Optimize.BFGS(neb)
                    case "fire":
                        optimizer = ASE_Optimize.FIRE(neb)
                    case "mdmin":
                        optimizer = ASE_Optimize.MDMin(neb)
                    case "gpmin":
                        optimizer = ASE_Optimize.GPMin(neb)
                    case _:
                        raise RuntimeError(
                            "Don't recognize NEB optimizer '" + P["neb optimizer"] + "'"
                        )
                self.optimizer = optimizer
                converged = optimizer.run(fmax=convergence, steps=max_steps)
        except Exception as err:  # noqa: F841
            self.logger.error("Caught exception: ", err)
            error = err
            raise_exception = True
        else:
            # Ran OK, so plot final rxn path and set the default structure to the best
            self.analyze(step=self._step)

            printer.normal("")
            if converged:
                text = f"The NEB calculation converged in {self._step} iterations. "
            else:
                text = (
                    f"The NEB calculation did not converge in {self._step} iterations. "
                )
            printer.normal(__(text, indent=indent))
        finally:
            toc = time.perf_counter_ns()
            self._results["t_elapsed"] = round((toc - tic) * 1.0e-9, 3)

            # Print the results
            # self.analyze()

            # Store results to db, variables, tables, and json as requested
            self.store_results(
                configuration=self._working_configuration,
                data=self._results,
            )
            self.optimizer = None

        # Clean up the subdirectories
        neb_dir = wd / "neb"
        if raise_exception:
            keep = P["on error"]
            if keep == "delete all subdirectories":
                shutil.rmtree(neb_dir)
            elif keep == "keep last subdirectory":
                subdirectories = neb_dir.glob("step_*")
                subdirectories = sorted(subdirectories)
                for subdirectory in subdirectories[:-1]:
                    shutil.rmtree(subdirectory)
            raise error from None
        else:
            keep = P["on success"]
            if keep == "delete all subdirectories":
                shutil.rmtree(neb_dir)
            elif keep == "keep last subdirectory":
                subdirectories = neb_dir.glob("step_*")
                subdirectories = sorted(subdirectories)
                for subdirectory in subdirectories[:-1]:
                    shutil.rmtree(subdirectory)

        if not converged and not P["continue if not converged"]:
            raise RuntimeError(
                f"The NEB calculation did not converge in {self._step} iterations"
            )

    def plot_path(self, step, plot_path, path, energies, fit_path, fit_energies, lines):
        """Plot the reaction path the to file 'plot'.

        Parameters
        ----------
        plot : pathlib.Path
            The file to write the plot to
        path : [float]
            The x points of the structures on the path
        energies : [float]
            The y points, or energies, of the structures
        fit_path : np.ndarray
            x points for the fit to the path
        fit_energies : np.ndarray
            y points for the fit path
        lines : [(float, float)]
            Tangent lines at each structure
        """

        figure = self.create_figure(
            module_path=(self.__module__.split(".")[0], "seamm"),
            template="line.graph_template",
            title=f"Reaction Path (step {step})",
        )

        plot = figure.add_plot("rxn_path")

        x_axis = plot.add_axis("x", label="Position along path")
        y_axis = plot.add_axis("y", label="E (kJ/mol)", anchor=x_axis)
        x_axis.anchor = y_axis

        plot.add_trace(
            x_axis=x_axis,
            y_axis=y_axis,
            name="Rxn Path",
            x=path,
            xlabel="position",
            y=energies,
            ylabel="E",
            yunits="kJ/mol",
            color="#4dbd74",
        )

        if fit_path is not None and fit_energies is not None:
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name="Fit",
                x=list(fit_path),
                xlabel="position",
                y=list(fit_energies),
                ylabel="Efit",
                yunits="kJ/mol",
                color="#000000",
            )

        figure.grid_plots("rxn_path")

        figure.dump(plot_path)

        write_html = "html" in self.options and self.options["html"]
        if write_html:
            figure.template = "line.html_template"
            figure.dump(plot_path.with_suffix(".html"))

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

        # Get the values of the parameters, dereferencing any variables
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Have to fix formatting for printing...
        PP = dict(P)
        for key in PP:
            if isinstance(PP[key], units_class):
                PP[key] = "{:~P}".format(PP[key])

        # Print what we are doing
        printer.important(__(self.description_text(P, short=True), indent=self.indent))

        # Create the directory
        wd = Path(self.directory)
        wd.mkdir(parents=True, exist_ok=True)

        # Do what we are asked
        approach = P["approach"].lower()
        if "interpolat" in approach:
            self.interpolate(P)
        elif "neb" in approach or "nudge" in approach:
            if "auto" in P["neb method"].lower():
                self.auto_neb(P)
            else:
                self.neb(P)
        else:
            raise RuntimeError("Don't understand approach '" + P["approach"] + "'")

        return next_node

    def set_id(self, node_id=()):
        """Sequentially number the subnodes"""
        self.logger.debug("Setting ids for subflowchart {}".format(self))
        if self.visited:
            return None
        else:
            self.visited = True
            self._id = node_id
            self.set_subids(self._id)
            return self.next()

    def set_subids(self, node_id=()):
        """Set the ids of the nodes in the subflowchart"""
        node = self.subflowchart.get_node("1").next()
        n = 1
        while node is not None:
            node = node.set_id((*node_id, str(n)))
            n += 1

    def write_structure(self, configuration, path):
        """Write the structure to disk for viewing.

        Parameters
        ----------
        configuration : molsystem._Configuration
            The structure to write.
        path : pathlib.Path or str
            The file to write to. The extension will be forced to mmcif and cif
        """
        if configuration.n_atoms > 0:
            # MMCIF file has bonds
            text = None
            try:
                text = configuration.to_mmcif_text()
            except Exception:
                message = "Error creating the mmcif file\n\n" + traceback.format_exc()
                self.logger.critical(message)

            if text is not None:
                Path(path).with_suffix(".mmcif").write_text(text)

            # CIF file has cell
            if configuration.periodicity == 3:
                text = None
                try:
                    text = configuration.to_cif_text()
                except Exception:
                    message = "Error creating the cif file\n\n" + traceback.format_exc()
                    self.logger.critical(message)

                if text is not None:
                    Path(path).with_suffix(".cif").write_text(text)

    def _write_sdf(self, path, configurations):
        """Write the configurations to an SDF file.

        Parameters
        ----------
        path : pathlib.Path or str
            The file to write
        configurations : molsystem._Configuration
            The configurations to write.
        """
        read_structure_step.write(
            str(path),
            configurations,
            extension=".sdf",
            remove_hydrogens=False,
            printer=None,
            references=self.references,
            bibliography=self._bibliography,
        )
