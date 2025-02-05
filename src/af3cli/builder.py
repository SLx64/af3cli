from typing import Self

from .input import InputFile
from .seqid import IDRegister
from .ligand import Ligand
from .sequence import Sequence
from .bond import Bond
from .sequence import Sequence, fasta2seq


class InputBuilder(object):
    """
    Constructs and configures an InputFile object.

    The InputBuilder class provides utility methods to set properties and add
    components to an InputFile object, enabling the creation of a fully
    configured input file for further processing. Each method modifies the
    InputFile instance held by the builder and returns the builder itself,
    allowing for method chaining.

    Attributes
    ----------
    _afinput : InputFile
        The input file being constructed and configured.
    _id_register : IDRegister
        The ID register object managing unique ID allocation.
    """
    def __init__(self, afinput: InputFile | None = None):
        if afinput is None:
            afinput = InputFile()
        self._afinput: InputFile = afinput
        self._id_register: IDRegister = IDRegister()

    def attach(self, afinput: InputFile) -> Self:
        """
        Attaches an `InputFile` object to the current instance.

        Parameters
        ----------
        afinput : InputFile
            The input file to be attached to the object.

        Returns
        -------
        self : object
            Returns the instance of the current object with the attached input
            file, allowing method chaining.
        """
        self._afinput = afinput
        return self

    def reset_ids(self) -> Self:
        """
        Resets all IDs of ligands and sequences within the current `Input` object.

        Returns
        -------
        Self
            Returns the instance itself to allow method chaining.
        """
        self._afinput.reset_ids()
        return self

    def build(self) -> InputFile:
        """
        Builds and retrieves the constructed `InputFile` instance.

        Returns
        -------
        InputFile
            The constructed and finalized `InputFile` instance.
        """
        return self._afinput

    def set_name(self, name: str) -> Self:
        """
        Sets the job name attribute in the `InputFile` instance.

        Parameters
        ----------
        name : str
            The new job name to set for the instance.

        Returns
        -------
        Self
            Returns the instance itself to allow method chaining.
        """
        self._afinput.name = name
        return self

    def set_seeds(self, seeds: list[int]) -> Self:
        """
        Sets the seeds attribute in the `InputFile` instance.

        Parameters
        ----------
        seeds : list of int
            A list of integers to be used as seeds.

        Returns
        -------
        Self
            Returns the instance of the object, allowing for method chaining.
        """
        self._afinput.seeds = set(seeds)
        return self

    def set_version(self, version: int) -> Self:
        """
        Sets the version attribute in the `InputFile` instance.

        Parameters
        ----------
        version : int
            The version number to set.

        Returns
        -------
        Self
            Returns the current instance of the object to allow method chaining.
        """
        self._afinput.version = version
        return self

    def set_dialect(self, dialect: str) -> Self:
        """
        Sets the dialect attribute in the `InputFile` instance.

        Parameters
        ----------
        dialect : str
            The dialect to set.

        Returns
        -------
        Self
            Returns the current instance of the object to allow method chaining.
        """
        self._afinput.dialect = dialect
        return self

    def add_ligand(self, ligand: Ligand) -> Self:
        """
        Adds a ligand to the collection of ligands in the `InputFile` instance.

        Parameters
        ----------
        ligand : Ligand
            The ligand to be added to the collection.

        Returns
        -------
        Self
            Returns the current instance of the object to allow method chaining.
        """
        self._afinput.ligands.append(ligand)
        return self

    def add_sequence(self, sequence: Sequence, 
                     seq_str: str = None, 
                     fasta_filename: str = None, 
                     num: int = 1 ) -> Self:
        """
        Adds a sequence to the list of sequences in the `InputFile` instance.

        Parameters
        ----------
        sequence : Sequence
            The sequence to be added to the list of sequences.
        sequence_str : str, optional
            The string representation of the sequence to be added.
        fasta_file : str, optional
            The path to the FASTA file containing the sequence to be added.
        num : int, optional
            The number of sequences to add. Defaults to 1.

        Returns
        -------
        Self
            Returns the current instance of the object to allow method chaining.
        """
        if sequence is None and seq_str is not None:
            sequence = Sequence(seq_str=seq_str, num=num)
        elif sequence is None and fasta_filename is not None:
            sequence =  next(fasta2seq(filename=fasta_filename, num=num))

        self._afinput.sequences.append(sequence)
        return self

    def add_bonded_atom_pair(self, bond: Bond) -> Self:
        """
        Add a bonded atom pair to the corresponding list in the `InputFile` instance.

        Parameters
        ----------
        bond : Bond
            An instance of the `Bond` class representing the bonded atom pair.

        Returns
        -------
        Self
            Returns the current instance of the object to allow method chaining.
        """
        self._afinput.bonded_atoms.append(bond)
        return self

    def set_user_ccd(self, user_ccd: str) -> Self:
        """
        Sets the userccd attribute as string in the `InputFile` instance.

        Parameters
        ----------
        user_ccd : str
            The CCD file content to be assigned to the userccd entry.

        Returns
        -------
        Self
            Returns the current instance of the object to allow method chaining.
        """
        self._afinput.user_ccd = user_ccd
        return self
