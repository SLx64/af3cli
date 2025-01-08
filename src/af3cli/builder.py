from typing import Self

from .input import InputFile
from .seqid import IDRegister
from .ligand import Ligand


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
    def __init__(self):
        self._afinput: InputFile = InputFile()
        self._id_register: IDRegister = IDRegister()

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
        self._afinput.seeds = seeds
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
