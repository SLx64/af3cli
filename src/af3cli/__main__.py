from __future__ import annotations

import pprint
import sys
import random
import logging
from typing import Callable, Self
from abc import ABCMeta, abstractmethod

import fire

from .builder import InputFile, InputBuilder
from .ligand import Ligand, LigandType, sdf2smiles
from .bond import Bond
from .sequence import Sequence, SequenceType
from .sequence import ProteinSequence, DNASequence, RNASequence
from .sequence import Template, TemplateType, MSA
from .sequence import (Modification, NucleotideModification,
                       ResidueModification)
from .sequence import read_fasta, fasta2seq
from .sequence import is_valid_sequence

# CONSTANTS
MAX_RANDOM_SEED: int = 99999
DEFAULT_FILENAME: str = "input.json"
CLI_NAME: str = "af3cli"

# LOGGER SETUP
logger = logging.getLogger("AF3 CLI")
logger.setLevel(logging.WARNING)

console_handler = logging.StreamHandler()
console_handler.setLevel(logging.WARNING)

formatter = logging.Formatter('%(name)s [%(levelname)s]: %(message)s')
console_handler.setFormatter(formatter)

logger.addHandler(console_handler)


def exit_on_error(msg: str) -> None:
    """
    Logs an error message and terminates the program.

    This function is used to report a critical error by logging the provided
    message using the logger and then terminating the program execution with
    a system exit.

    Parameters
    ----------
    msg : str
        The error message to log before exiting the program.

    Returns
    -------
    None
        This function does not return; instead, it terminates the program.
    """
    logger.error(msg)
    sys.exit(1)


def random_int_list(start: int, end: int, num: int) -> list[int]:
    """
    Generate a random list of unique integers within a specified range.

    This function generates a list of unique integers chosen randomly from
    a specified range `[start, end]` (inclusive). The total number of integers
    to be chosen is specified by the `num` parameter. The function ensures
    that the generated list contains no duplicate values.

    Parameters
    ----------
    start : int
        The inclusive lower bound of the range.
    end : int
        The inclusive upper bound of the range. Must be greater than `start`.
    num : int
        The number of unique integers to generate. Must not exceed
        the size of the range (`end - start + 1`).

    Returns
    -------
    list of int
        A list of unique integers randomly generated within the specified range.
    """
    if end <= start or num > end-start:
        raise ValueError("Invalid range specified.")
    result = set()
    while len(result) < num:
        result.add(random.randint(start, end))
    return list(result)


def read_sdf_file(filename: str) -> list[str]:
    """
    Reads a Structure Data File (SDF) and converts the molecules into SMILES format.

    Parameters
    ----------
    filename : str
        The path to the SDF file that needs to be read.

    Returns
    -------
    list of str
        A list of SMILES strings representing the molecules contained in the
        specified SDF file.

    Raises
    ------
    SystemExit
        If RDKit is not installed on the system or file is not found.
    """
    try:
        smiles = []
        smi_generator = sdf2smiles(filename)
        for smi in smi_generator:
            if smi is None:
                logger.warning("Failed to read molecule from SDF file.")
                continue
            smiles.append(smi)
        logger.info(f"Read {len(smiles)} molecules from SDF file.")
        return smiles
    except FileNotFoundError as e:
        exit_on_error(f"SDF file not found: {e}")
    except ImportError as e:
        exit_on_error(e.msg)
    except Exception as e:
        exit_on_error(f"Failed to read SDF file: {e}")


def read_fasta_entry(filename: str) -> str:
    """
    Reads a single entry from a FASTA file and returns its sequence string.

    Parameters
    ----------
    filename : str
        The path to the FASTA file to be read.

    Returns
    -------
    str
        The sequence string from the first valid entry in the FASTA file.

    Raises
    ------
    SystemExit
        If Biopython is not installed or if the specified FASTA file
        does not exist / contains no valid sequences.
    """
    try:
        fasta_file = read_fasta(filename)
        seq_name, seq_str = next(fasta_file)
        if seq_str is None:
            exit_on_error("No valid sequence found in FASTA file.")
        return seq_str
    except FileNotFoundError as e:
        exit_on_error(f"FASTA file not found: {e}")
    except StopIteration as e:
        exit_on_error(f"No sequence found in FASTA file. {e}")
    except ImportError as e:
        exit_on_error(e.msg)
    except Exception as e:
        exit_on_error(f"Failed to read FASTA file: {e}")


def read_fasta_file(filename: str) -> list[Sequence]:
    """
    Reads sequences from a FASTA file and returns a list of sequences.

    Parameters
    ----------
    filename : str
        Path to the FASTA file to read.

    Returns
    -------
    list of Sequence
        A list containing the parsed sequences from the specified FASTA file.

    Raises
    ------
    SystemExit
        If Biopython is not installed or if the specified FASTA file
        does not exist.
    """
    try:
        fasta_file = fasta2seq(filename)
        sequences = []
        for entry in fasta_file:
            if entry is None:
                logger.warning("Failed to read sequence from FASTA file.")
                continue
            sequences.append(entry)
        if len(sequences) == 0:
            logger.warning("No valid sequences found in FASTA file.")
        return sequences
    except FileNotFoundError as e:
        exit_on_error(f"FASTA file not found: {e}")
    except ImportError as e:
        exit_on_error(e.msg)
    except Exception as e:
        exit_on_error(f"Failed to read FASTA file: {e}")


def read_file_to_str(filename) -> str:
    """
    Reads the content of a specified file and returns it as a string.

    In the AlphaFold3 input file, some variables can be specified as
    string or file path. Since adding large file content as a command line
    argument is not practical, this function is used to read the content of
    the corresponding file if necessary.

    Parameters
    ----------
    filename : str
        The path to the file to be read.

    Returns
    -------
    str
        The content of the file as a string.
    """
    logger.info(f"Reading file '{filename}'")
    with open(filename, "r") as f:
        return f.read()


def ensure_int_list(value: int | tuple | list | None) -> list[int]:
    """
    Converts a given input to a list of integers. The function ensures the
    input is consistently returned as a list of integers, regardless of
    whether the input is a single integer, a tuple, or None. If None is
    passed, an empty list is returned.

    Parameters
    ----------
    value : int, tuple, or None
        The input value to be converted to a list of integers.

    Returns
    -------
    list of int
        A list of integers derived from the input value. If the input is `None`,
        an empty list is returned.
    """
    if value is None:
        return []
    elif isinstance(value, tuple):
        return list(value)
    elif isinstance(value, int):
        return [value]
    return value


def ensure_opt_str_list(value: str | tuple | list | None) -> list[str] | None:
    """
    Ensures the provided input is converted to a list of strings or None, since
     --ids="A,B" (tuple) and --ids="[A,B]" (list) are both valid inputs.

    Parameters
    ----------
    value : str, tuple, list, or None
        The input value to be processed. This can be of type string, tuple, list, or
        None.

    Returns
    -------
    list, None, or input value
        Returns the input as a list of strings if it is a string or tuple. If the
        input is `None`, it returns `None`. All other inputs are returned unchanged.

    """
    if value is None:
        return None
    elif isinstance(value, tuple):
        return list(value)
    elif isinstance(value, str):
        return [value]
    return value


def hide_from_cli(function: Callable) -> Callable:
    """
    A decorator to mark a function as hidden from a CLI interface to prevent
    an automatic inclusion as command line argument or subcommand.

    Parameters
    ----------
    function : callable
        The function to decorate and mark as hidden from the CLI interface.

    Returns
    -------
    callable
        The same function passed as input, with the `_hidden` attribute set
        to `True`.
    """
    function._hidden = True
    return function


class CommandBase(object, metaclass=ABCMeta):
    """
    Base class for command-like objects.

    This class serves as a foundational structure for implementing command-like objects
    It includes mechanisms to hide certain attributes from being accessed via the
    command-line interface (CLI) and an abstract method `finalize` to enforce implementation
    in derived classes.
    """
    def __dir__(self) -> list[str]:
        """
        Filters the attributes of an object to exclude attributes marked as hidden.

        This method overrides the default `__dir__` implementation to provide a
        filtered list of visible attributes for the object. Attributes with `_hidden`
        set to `True` will be excluded from the resulting list.

        Returns
        -------
        list of str
            A list of attribute names that are not marked as hidden.

        """
        attributes = super().__dir__()
        return [attr for attr in attributes
                if not getattr(getattr(self, attr, None), "_hidden", False)]

    def __call__(self) -> None:
        """
        This method is invoked when a subcommand completes its execution.
        """
        self.finalize()

    @abstractmethod
    def finalize(self) -> None:
        """
        This abstract method is required to be implemented by any subclass. It is
        intended to encapsulate the finalization logic necessary to properly wrap
        up the processing or operations associated for the corresponding command.

        The function will be automatically invoked when the command is executed.
        """
        pass


class CLICommand(CommandBase, metaclass=ABCMeta):
    """
    Base class for command-like objects.

    This class serves as a foundational structure for implementing command-like objects
    It includes mechanisms to hide certain attributes from being accessed via the
    command-line interface (CLI)
    and an abstract method `finalize` to enforce implementation in derived classes.

    Attributes
    ----------
    _parent : CLICommand or None
        Stores the parent command object or None if no parent is set.
    """
    def __init__(self):
        self._parent: CLI | None = None

    @hide_from_cli
    def get_parent(self) -> CLI | None:
        """
        Gets the parent command object.

        This method retrieves the parent command object associated
        with the current instance. The parent is typically used
        to access the builder object in the top level CLI object.

        Returns
        -------
        CommandBase
            The parent command object.
        """
        return self._parent

    @hide_from_cli
    def set_parent(self, parent: CLI) -> Self:
        """
        Sets the parent command and returns the current instance.

        Parameters
        ----------
        parent : CommandBase
            The parent command to be assigned to the current instance.

        Returns
        -------
        CommandBase
            Returns the current command instance after assigning the parent.
        """
        self._parent = parent
        return self

    def __call__(self) -> CLI | None:
        """
        This method is invoked when a subcommand completes its execution.

        It ensures any necessary finalization tasks are performed by
        calling `finalize()` and then returns control back to the parent
        command by returning the parent object.
        """
        self.finalize()
        return self._parent


class SequenceCommand(CLICommand, metaclass=ABCMeta):
    """
    Represents a sequence command used for managing sequences in the AlphaFold3
    input file. This class allows for the addition of sequences, their
    modifications, multiple sequence alignment (MSA), and templates. It is
    intended to serve as a base class for more specific sequence-related commands.

    Attributes
    ----------
    _parent : CLI or None
        The parent object to which this sequence command belongs.
    _sequence_str : str or None
        The sequence string.
    _sequence_ids : list[str], str or None
        A list or single identifier(s) associated with the sequence.
    _sequence_num : int
        The number of sequences to be added. Will be overwritten if
        `_sequence_ids` is not None.
    _modifications : list[Modification] or None
        A list of modifications applied to the sequence.
    _templates : list[Template] or None
        A list of templates associated with the sequence.
    _msa : MSA or None
        The multiple sequence alignment (MSA) object associated with the sequence.
    """
    def __init__(self):
        super().__init__()
        self._parent: CLI | None = None
        self._sequence_str: str | None = None
        self._sequence_ids: list[str] | str | None = None
        self._sequence_num: int = 1
        self._modifications: list[Modification] | None = None
        self._templates: list[Template] | None = None
        self._msa: MSA | None = None

    def add(
        self,
        sequence: str,
        num: int = 1,
        ids: list[str] | None = None,
        fasta: bool = False
    ) -> Self:
        """
        Adds a sequence and optionally specifies either a number or a list of IDs for the
        sequence. These parameters should not be used simultaneously.

        Parameters
        ----------
        sequence : str
            A string representing the specific sequence to be added or a path
            to a FASTA file containing the sequence if `fasta` is set to True.

        num : int
            An integer that defines a number of sequences. Either
            `num` or `ids` can be specified, but not both.

        ids : list of str
            A list of strings representing unique identifiers for the sequence.
            Either `ids` or `num` can be specified, but not both.

        fasta: bool
            If `True`, the sequence will be interpreted as a FASTA file.

        Returns
        -------
        SequenceCommand
            The current instance of the object updated with the provided sequence and
            parameters.

        """
        ids = ensure_opt_str_list(ids)

        if fasta:
            self._sequence_str = read_fasta_entry(sequence)
        else:
            self._sequence_str = sequence

        if not is_valid_sequence(self.sequence_type(), self._sequence_str):
            exit_on_error(f"Invalid sequence for '{self.sequence_type().name}'")

        self._sequence_ids = ids
        self._sequence_num = num
        return self

    def modification(self, mod: str, pos: int) -> Self:
        """
        Subcommand for adding a specific modification at a given position.
        Supports different types of sequences, such as protein and nucleotide,
        appending the modification instance accordingly.

        Parameters
        ----------
        mod : str
            The modification string representing the CCD code.

        pos : int
            The position in the sequence where the modification is applied. It
            should be a 1-based index.

        Returns
        -------
        SequenceCommand
            Returns the current sequence command object after appending the
            modification, enabling method chaining.
        """
        if self._modifications is None:
            self._modifications = []

        if self.sequence_type() == SequenceType.PROTEIN:
            self._modifications.append(
                ResidueModification(mod_str=mod, mod_pos=pos)
            )
        else:
            self._modifications.append(
                NucleotideModification(mod_str=mod, mod_pos=pos)
            )
        return self

    def msa(
        self,
        paired: str | None = None,
        unpaired: str | None = None,
        pairedpath: str | None = None,
        unpairedpath: str | None = None
    ) -> Self:
        """
        Subcommand to add multiple sequence alignment (MSA).

        The string and path variants for paired and unpaired sequences are mutually exclusive.

        Parameters
        ----------
        paired : str
            A string representation of paired sequences to be aligned. This parameter is
            mutually exclusive with `pairedpath`.
        unpaired : str
            A string representation of unpaired sequences to be aligned. This parameter is
            mutually exclusive with `unpairedpath`.
        pairedpath : str
            The file path to paired sequences to be aligned. This parameter is mutually
            exclusive with `paired`.
        unpairedpath : str
            The file path to unpaired sequences to be aligned. This parameter is mutually
            exclusive with `unpaired`.

        Returns
        -------
        SequenceCommand
            An instance of the sequence command constructed with the MSA configuration based
            on the provided arguments.

        Raises
        ------
        SystemExit
            If both string and path variants for paired or unpaired sequences are provided,
            an error is raised, indicating they are mutually exclusive.
        """
        if ((paired is not None and pairedpath is not None) or
                (unpaired is not None and unpairedpath is not None)):
            exit_on_error("String and path variants "
                          "for MSA are mutually exclusive.")

        paired = paired or unpairedpath
        unpaired = unpaired or unpairedpath
        self._msa = MSA(
            paired=paired,
            unpaired=unpaired,
            paired_is_path=pairedpath is not None,
            unpaired_is_path=unpairedpath is not None
        )
        return self

    @abstractmethod
    def sequence_type(self) -> SequenceType:
        """
        Abstract method to retrieve the sequence type.

        This method must be implemented by all derived classes, providing the
        specific type of sequence that the implementation represents.

        Returns
        -------
        SequenceType
            Type of the sequence associated with the implementation.
        """
        pass


class ProteinCommand(SequenceCommand):
    """
    Represents a command to work with protein sequences.

    This class extends the functionality of `SequenceCommand` by adding
    specialized operations for handling protein-specific data, including
    templates, modifications, and multiple sequence alignment (MSA).

    """
    def __init__(self):
        super().__init__()

        self._templates: list[Template] = []

    def template(
        self,
        mmcif: str,
        qidx: list[int] | None = None,
        tidx: list[int] | None = None,
        read: bool = False
    ) -> Self:
        """
        Subcommand to add a template from an mmCIF file.

        This function provides the ability to add templates based on the mmCIF format,
        as file path only. The file content can be read into a string if the parameter
        `read` is set to True.

        Parameters
        ----------
        mmcif : str
            The file path to the mmCIF file.
        qidx : list of int
            A list of query indices. Defaults to None.
        tidx : list of int
            A list of template indices. Defaults to None.
        read : bool
            If True, the `mmcif` file content is read into a string. Defaults to False.

        Returns
        -------
        ProteinCommand
            The current instance of the class, allowing for method chaining.
        """
        if read:
            mmcif = read_file_to_str(mmcif)
            template_type = TemplateType.STRING
        else:
            template_type = TemplateType.FILE

        self._templates.append(
            Template(
                mmcif=mmcif,
                template_type=template_type,
                qidx=ensure_int_list(qidx),
                tidx=ensure_int_list(tidx)
            )
        )
        return self

    @hide_from_cli
    def sequence_type(self) -> SequenceType:
        """
        Sets the sequence type to protein.

        Returns
        -------
        SequenceType
            Returns the sequence type as a `SequenceType.PROTEIN` enum value.

        """
        return SequenceType.PROTEIN

    @hide_from_cli
    def finalize(self) -> None:
        """
        Finalize and add the constructed sequence to the parent builder.
        """
        sequence = ProteinSequence(
            seq_str=self._sequence_str,
            num=self._sequence_num,
            seq_id=self._sequence_ids,
            modifications=self._modifications or None,
            templates=self._templates or None,
            msa=self._msa
        )
        self._parent.builder().add_sequence(sequence)
        logger.info(f"Adding sequence {sequence}")

        self._templates = []
        self._modifications = None
        self._msa = None


class DNACommand(SequenceCommand):
    """
    Represents a command to work with DNA sequences.

    This class extends the functionality of `SequenceCommand` by adding
    specialized operations for handling DNA-specific data, including
    templates, modifications, and multiple sequence alignment (MSA).
    """
    def __init__(self):
        super().__init__()
        self._rev_complement: bool = False

    def add(
        self,
        sequence: str,
        num: int = 1,
        ids: list[str] | None = None,
        fasta: bool = False,
        complement: bool = False
    ) -> Self:
        """
        Adds a given sequence to the current object with options for specifying
        the number of times to add, associated IDs, and whether to use FASTA
        format or complementary sequences.

        Parameters
        ----------
        sequence : str
            The sequence to be added.
        num : int
            The number of times the sequence should be added. Default is 1.
        ids : list[str] or None, optional
            A list of string identifiers associated with the sequence. Default
            is None.
        fasta : bool, optional
            A flag to determine if the sequence should be added in FASTA format.
            Default is False.
        complement : bool, optional
            A flag to indicate if the reverse complement of the sequence should
            also be added. Modifications and IDs will be ignored. Default is False.

        Returns
        -------
        Self
            The updated instance of the calling object.
        """
        self._rev_complement = complement
        return super().add(sequence, num, ids, fasta)

    @hide_from_cli
    def sequence_type(self) -> SequenceType:
        """
        Sets the sequence type to DNA.

        Returns
        -------
        SequenceType
            Returns the sequence type as a `SequenceType.DNA` enum value.

        """
        return SequenceType.DNA

    @hide_from_cli
    def msa(self):
        raise NotImplementedError("MSA is not available for DNA sequences.")

    @hide_from_cli
    def finalize(self) -> None:
        """
        Finalize and add the constructed sequence to the parent builder.
        """
        sequence = DNASequence(
            seq_str=self._sequence_str,
            num=self._sequence_num,
            seq_id=self._sequence_ids,
            modifications=self._modifications or None
        )
        self._parent.builder().add_sequence(sequence)
        logger.info(f"Adding sequence {sequence}")

        if self._rev_complement:
            complement = sequence.reverse_complement()
            self._parent.builder().add_sequence(complement)
            logger.info(f"Adding sequence {complement}")
            if self._sequence_ids is not None:
                logger.warning("ID assignment for complementary sequences "
                               "is not supported. Please add the complementary "
                               "sequence separately.")
            if self._modifications is not None:
                logger.warning("Modification assignment for complementary "
                               "sequences is not supported. Please add the"
                               "complementary sequence separately.")

        self._modifications = None
        self._rev_complement = False


class RNACommand(SequenceCommand):
    """
    Represents a command to work with RNA sequences.

    This class extends the functionality of `SequenceCommand` by adding
    specialized operations for handling RNA-specific data, including
    templates, modifications, and multiple sequence alignment (MSA).
    """
    def __init__(self):
        super().__init__()

    @hide_from_cli
    def sequence_type(self) -> SequenceType:
        """
        Sets the sequence type to RNA.

        Returns
        -------
        SequenceType
            Returns the sequence type as a `SequenceType.RNA` enum value.
         """
        return SequenceType.RNA

    @hide_from_cli
    def finalize(self) -> None:
        """
        Finalize and add the constructed sequence to the parent builder.
        """
        sequence = RNASequence(
            seq_str=self._sequence_str,
            num=self._sequence_num,
            seq_id=self._sequence_ids,
            modifications=self._modifications or None,
            msa=self._msa
        )
        self._parent.builder().add_sequence(sequence)
        logger.info(f"Adding sequence {sequence}")

        self._modifications = None
        self._msa = None


class LigandCommand(CLICommand):
    """
    Represents a command to add ligands.

    This class allows the user to define ligands using different input types such as
    SMILES strings, CCD codes, or SDF files. It ensures that only one input type is
    used per command. The ligands are added into a parent object's builder.

    Attributes
    ----------
    _ligands : list[Ligand] or None
        A list of Ligand objects representing the added ligands.
    """
    def __init__(self):
        super().__init__()
        self._ligands: list[Ligand] | None = None

    def add(
        self,
        smiles: str | None = None,
        ccd: list[str] | str | None = None,
        sdf: str | None = None,
        num: int = 1,
        ids: list[str] | None = None
    ) -> Self:
        """
        Adds a ligand entry to the command with proper validation and formatting.

        This method allows adding ligands to a command instance using one of the
        supported formats (SMILES, CCD, or SDF). It ensures that only one ligand
        variant is provided at a time and handles cases where either a number or
        specific IDs are supplied for the ligands. Additionally, in the context of
        SDF files, it validates that the number of IDs corresponds correctly to the
        entries within the file.

        Parameters
        ----------
        smiles : str
            SMILES representation of the ligand if provided. Must not be used
            simultaneously with `ccd` or `sdf`.

        ccd : list[str] or str
            CCD identifier of the ligand if provided. Must not be used
            simultaneously with `smiles` or `sdf`.

        sdf : str
            Path to the SDF file containing the ligand data. Must not be used
            simultaneously with `smiles` or `ccd`.

        num : int
            Number of ligands to process. Will be overwritten if `ids` is
            specified. Defaults to 1.

        ids : list[str]
            List of specific IDs corresponding to the ligands. Used only with
            `smiles`, `ccd`, or entries in an SDF file.

        Returns
        -------
        LigandCommand
            The instance of the calling class with the updated ligand list.
        """
        if self._ligands is not None:
            exit_on_error("Command 'add' was called more than once.")

        if sum(arg is not None for arg in [smiles, ccd, sdf]) != 1:
            exit_on_error("Only one variant should be provided for ligands")

        # create an empty list of ligands since multiple ligands can be
        # added at once with an SDF file.
        self._ligands = []
        ids = ensure_opt_str_list(ids)

        if sdf is not None:
            sdf_smiles = read_sdf_file(sdf)

            # The number of ligands in the SDF file should match
            # the correponding id values if specified, otherwise the
            # ids cannot be assigned
            if ids is not None:
                logger.info("Explicit ID assignment is not needed for SDF "
                            "as each entry is converted to SMILES.")
                if len(ids) != len(sdf_smiles):
                    exit_on_error("Number of ids does not "
                                  "match number of SDF entries.")
            else:
                # ids can be automatically assigned
                ids = [None for _ in range(len(sdf_smiles))]

            for k, smi in zip(ids, sdf_smiles):
                k = [k] if k is not None else None
                self._ligands.append(
                    Ligand(
                        ligand_type=LigandType.SMILES,
                        ligand_value=smi,
                        num=num,
                        seq_id=k
                    )
                )
            return self

        if isinstance(smiles, list | tuple):
            exit_on_error("Only a single SMILES string can be provided.")

        ligand_str = smiles or ensure_opt_str_list(ccd)
        ligand_type = LigandType.SMILES if smiles else LigandType.CCD

        # append the single ligand to a list to unify the finalization code
        self._ligands.append(
            Ligand(
                ligand_type=ligand_type,
                ligand_value=ligand_str,
                num=num,
                seq_id=ids
            )
        )
        return self

    @hide_from_cli
    def finalize(self) -> None:
        """
        Finalizes the ligand addition process in the builder context.

        This method iterates over all ligands stored in the instance and
        adds them to the associated parent's builder. Once the process is complete,
        it resets the `_ligands` attribute to `None` since multiple ligand commands
        are allowed.
        """
        if self._ligands is None:
            exit_on_error("No ligands were added.")

        for ligand in self._ligands:
            logger.info(f"Adding ligand {ligand}")
            self._parent.builder().add_ligand(ligand)
        self._ligands = None


class CLI(CommandBase):
    """
    CLI to manage AlphaFold3 input configurations.

    This CLI provides various operations to manage and create input configurations for
    AlphaFold jobs, including handling commands for proteins, DNA, RNA, and ligands. It
    provides utilities to manage the basic configurations and random seed setups.
    """
    def __init__(self):
        super().__init__()
        self._builder: InputBuilder = InputBuilder()
        self._filename: str = DEFAULT_FILENAME

        self._debug_print: bool = False

        self.protein: ProteinCommand = ProteinCommand().set_parent(self)
        self.dna: DNACommand = DNACommand().set_parent(self)
        self.rna: RNACommand = RNACommand().set_parent(self)
        self.ligand: LigandCommand = LigandCommand().set_parent(self)

    def debug(self, verbose: bool = False, show: bool = False) -> Self:
        """
        Command to set the verbosity and debug output level for logging.

        This method adjusts the logging levels for the logger and console handler based on the
        `verbose` flag. If `verbose` is `True`, the logging level is set to `DEBUG`; otherwise,
        it defaults to `WARNING`. Additionally, the `show` parameter determines whether the
        AlphaFold3 input file is printed to the console or written to a file.

        Parameters
        ----------
        verbose : bool
            Determines whether detailed debug logging is enabled. If set to True, the logging
            level is set to DEBUG; if False, it is set to WARNING. Defaults to False.

        show : bool
            Enables or disables debug printing of the final AlphaFold3 input file.
            When set to True, debug printing is turned on. Defaults to False.

        Returns
        -------
        CLI
            Returns the same instance of the class to enable method chaining.
        """
        self._debug_print = show

        if verbose:
            logger.setLevel(logging.DEBUG)
            console_handler.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.WARNING)
            console_handler.setLevel(logging.WARNING)
        return self

    def config(
        self,
        filename: str = DEFAULT_FILENAME,
        jobname: str = "job",
        version: int = 1,
        dialect: str = "alphafold3"
    ) -> Self:
        """
        Command to add basic information to the AlphaFold3 input file,
        like the job name, version, and dialect, as well as the file name.

        Parameters
        ----------
        filename : str, default="input.json"
            The name of the output file.
        jobname : str, default="job"
            The name of the job for which the configuration is set.
        version : int, default=1
            The version number of the AlphaFold3 configuration.
        dialect : str, default="alphafold3"
            The AlphaFold3 dialect specification.

        Returns
        -------
        CLI
            Returns the same instance of the class to enable method chaining.
        """
        self._filename = filename
        self._builder.set_name(jobname)
        self._builder.set_version(version)
        self._builder.set_dialect(dialect)
        return self

    def merge(
        self,
        filename: str,
        noreset: bool = False,
        userccd: bool = False,
        bonds: bool = False,
        seeds: bool = False
    ) -> CLI:
        """
        Merges the current input object with a specified input file.

        Parameters
        ----------
        filename : str
            The path to the input JSON file to be merged with the current input.
        noreset : bool, optional
            If True, skips resetting IDs while merging inputs
            (CAUTION: this might result in clashes with existing IDs).
            Defaults to False.
        userccd : bool, optional
            If True, overrides user-defined CCD while merging.
            Defaults to False.
        bonds : bool, optional
            If True, merges bonded atoms data.
            Defaults to False.
        seeds : bool, optional
            If True, merges seeds and removes duplicates.
            Defaults to False.

        Returns
        -------
        CLI
            Returns the same instance of the class to enable method chaining.
        """
        try:
            other_input = InputFile.read(filename)
            if noreset:
                logger.warning("Skipping reset might cause ID clashes.")

            curr_input = self._builder.build()
            curr_input.merge(
                other_input,
                reset=not noreset,
                userccd=userccd,
                bonded_atoms=bonds,
                seeds=seeds
            )
        except Exception as e:
            exit_on_error(f"Failed to process existing input file: {filename}\n{e}")
        return self

    def fasta(self, filename: str) -> Self:
        """
        Command to process and incorporate sequences from a FASTA file.

        Parameters
        ----------
        filename : str
            The file path of the FASTA file to be read and processed.

        Returns
        -------
        CLI
            Returns the same instance of the class to enable method chaining.
        """
        for seq in read_fasta_file(filename):
            self._builder.add_sequence(seq)
        return self

    def ccd(self, filename: str) -> Self:
        """
        Command to process and incorporate user-provided CCD data.

        This method reads the CCD data from a specified file and
        adds it to the input configuration.

        Parameters
        ----------
        filename : str
            The file path of the CCD file to be read and processed.

        Returns
        -------
        CLI
            Returns the same instance of the class to enable method chaining.
        """
        user_ccd = read_file_to_str(filename)
        self._builder.set_user_ccd(user_ccd)
        return self

    def seeds(self, values: list[int] = None, num: int = None) -> Self:
        """
        Command to add new seed values to the configuration.

        This method allows the user to either specify a list of seed values directly
        or a number to generate a random list of seeds.

        Parameters
        ----------
        values : list[int]
            A list of integer seed values. If provided, overrides the `num` parameter.
            Must be explicitly specified if `num` is not provided.
        num : int
            The number of random integer seed values to be generated. This parameter
            is ignored if `values` is provided.

        Returns
        -------
        CLI
            Returns the same instance of the class to enable method chaining.

        Notes
        -----
        - The maximum random seed value used during random generation is determined
          by the `MAX_RANDOM_SEED` constant within the program scope.
        """
        if sum(arg is not None for arg in [values, num]) != 1:
            exit_on_error("Specify either 'values' or 'num'.")
        if num is not None:
            values = random_int_list(1, MAX_RANDOM_SEED, num=num)

        # convert to a list to conform the JSON format
        values = ensure_int_list(values)
        if len(values):
            self._builder.set_seeds(values)
        return self

    def bond(self, add: str) -> Self:
        """
        Command to add a new bond between two atoms.

        This method takes a string representation, converts it into a
        Bond object, and adds the bonded atom pair to the current
        input file using the builder. The sequence id, residue number and
        atom name are separated by a colon and both atoms are separated by
        a dash. Example: "1:A:CA-1:A:CB"

        Parameters
        ----------
        add : str
            A string representation of the bond to be added, which specifies
            the atom pair to be bonded.

        Returns
        -------
        CLI
            Returns the same instance of the class to enable method chaining.
        """
        b = Bond.from_string(add)
        self._builder.add_bonded_atom_pair(b)
        return self

    @hide_from_cli
    def builder(self) -> InputBuilder:
        """
        Returns the builder instance associated with this object.

        This method provides access to the builder object which can be used
        to construct and configure the necessary inputs from subcommands.

        Returns
        -------
        InputBuilder
            The builder instance used for constructing input files.
        """
        return self._builder

    @hide_from_cli
    def finalize(self) -> None:
        """
        Finalizes the construction and writing of the AF3 input file.

        This method builds the input file  and writes the file to the specified
        filename. If `_debug_print` is enabled, the generated file data is
        printed in a pretty-printed dictionary format. Otherwise, the method
        writes the file content in JSON format to the specified location.
        """
        af_input_file = self._builder.build()
        if self._debug_print:
            pp = pprint.PrettyPrinter(indent=4)
            pp.pprint(af_input_file.to_dict())
        else:
            af_input_file.write(self._filename)
            logger.info(f"Writing AF3 input file to '{self._filename}'")


def main() -> None:
    # enable printing of the help to the console
    fire.core.Display = lambda lines, out: print(*lines, file=out)

    if len(sys.argv) == 1:
        fire.Fire(CLI(), name=CLI_NAME, command=["--", "--help"])
    else:
        fire.Fire(CLI(), name=CLI_NAME)


if __name__ == '__main__':
    main()
