from __future__ import annotations

from copy import deepcopy

from .mixin import DictMixin
from .ligand import Ligand
from .bond import Bond
from .sequence import Sequence
from .seqid import IDRegister


class InputFile(DictMixin):
    """
    Represents an input file configuration for an AlphaFold3 job, including
    metadata, sequences and ligands. This class facilitates the preparation
    and conversion of input data to dictionary format for further
    processing or serialization as JSON file.

    Attributes
    ----------
    name : str
        Name of the job.
    version : int
        Version of the input file format.
    dialect : str
        Input file dialect.
    seeds : list of int
        List of seed integers.
    user_ccd : str or None
        Custom CCD structural data provided by the user, if any.
    """
    def __init__(
        self,
        name: str = "job",
        version: int = 1,
        dialect: str = "alphafold3",
        seeds: list[int] | None = None,
        user_ccd: str | None = None,
    ):
        self.name: str = name
        self.version: int = version
        self.dialect: str = dialect
        self.user_ccd: str = user_ccd
        self.seeds: set[int] = set()

        if seeds is None:
            seeds = [1]
        self.seeds.update(seeds)

        self.ligands: list[Ligand] = []
        self.bonded_atoms: list[Bond] = []
        self.sequences: list[Sequence] = []

        self._id_register: IDRegister = IDRegister()

    def _register_ids(self) -> None:
        """
        Registers unique IDs for sequences and ligands in the internal ID register.

        This method iterates through the objects stored in sequences and ligands,
        checking if each entry is already registered. If an entry is not yet
        registered, it retrieves the associated IDs by calling its `get_id`
        method. Each ID from the returned list is then registered into the
        internal `IDRegister`.

        Raises
        ------
        AttributeError
            If an item in sequences or ligands does not have `is_registered` or
            `get_id` attributes.
        """
        self._id_register.reset()
        self.clear_temporary_ids()
        for seqtype in [self.sequences, self.ligands]:
            for entry in seqtype:
                seq_ids = entry.get_id()
                if seq_ids is not None:
                    for seq_id in seq_ids:
                        self._id_register.register(seq_id)

    def _assign_ids(self) -> None:
        """
        Assign unique IDs to sequences and ligands that do not already have an ID.

        This method iterates over sequences and ligands, checking if they have an ID.
        If IDs are missing, it generates unique IDs for entries within each type and
        assigns these IDs.
        """
        for seqtype in [self.sequences, self.ligands]:
            for entry in seqtype:
                num_ids = entry.required_tmp_id_count()
                seq_ids = [self._id_register.generate() for _ in range(num_ids)]
                entry.set_temporary_id(seq_ids)

    def _prepare(self) -> None:
        self._register_ids()
        self._assign_ids()

    def reset_all_ids(self) -> None:
        """
        Resets the IDs of all entries in the object's sequences and ligands.
        """
        for seqtype in [self.sequences, self.ligands]:
            for entry in seqtype:
                entry.remove_id()
                entry.clear_temporary_id()

    def clear_temporary_ids(self) -> None:
        for seqtype in [self.sequences, self.ligands]:
            for entry in seqtype:
                entry.clear_temporary_id()

    def merge(
        self,
        other: InputFile,
        reset: bool = True,
        seeds: bool = False,
        bonded_atoms: bool = False,
        userccd: bool = False
    ) -> None:
        """
        Merges the content of another InputFile instance into the current
        instance.

        Notes
        -----
        No checks are performed to ensure that the IDs are handled correctly.
        Please use with caution when no reset is performed or when keeping
        bonded atoms.

        Parameters
        ----------
        other : InputFile
            An InputFile instance whose content will be merged into
            the current instance.
        reset : bool
            If True (default), resets the IDs of the content being merged
            to ensure unique identifiers.
        seeds : bool
            If True (default: False), merges the seeds from the `other`
            InputFile instance into the current instance.
        bonded_atoms : bool
            If True (default: False), appends bonded atom information from
            the `other` InputFile instance to the current instance.
        userccd : bool
            If True (default: False), overwrites the `user_ccd` attribute
            of the current instance with the value from the `other` InputFile
            instance.

        Returns
        -------
        None
            This method modifies the current instance in place; it does
            not return any value.
        """
        tmp_input = deepcopy(other)

        if reset:
            tmp_input.reset_all_ids()
        if seeds:
            self.seeds.update(tmp_input.seeds)
        if bonded_atoms:
            for bond in tmp_input.bonded_atoms:
                self.bonded_atoms.append(bond)
        if userccd:
            self.user_ccd = tmp_input.user_ccd

        for seq in tmp_input.sequences:
            self.sequences.append(seq)

        for lig in tmp_input.ligands:
            self.ligands.append(lig)

    def to_dict(self) -> dict:
        """
        Converts the object and its associated attributes to a dictionary representation.
        The created dictionary contains all necessary information for the AlphaFold3 input
        file. An ID will be assigned to all sequences and ligands if they do not already
        have one. This might result in an error if new sequences with duplicate IDs are
        added to the input file.

        Returns
        -------
        dict
            A dictionary containing all relevant attributes for running AlphaFold3 jobs.
        """
        self._prepare()

        content = dict()
        content["name"] = self.name
        content["version"] = self.version
        content["dialect"] = self.dialect
        content["modelSeeds"] = list(self.seeds)
        content["sequences"] = []

        for seqtype in [self.sequences, self.ligands]:
            for entry in seqtype:
                entry_dict = entry.to_dict()

                if self.version < 4:
                    for inner in entry_dict.values():
                        inner.pop("description", None)

                content["sequences"].append(entry_dict)

        if len(self.bonded_atoms):
            content["bondedAtomPairs"] = []
            for entry in self.bonded_atoms:
                content["bondedAtomPairs"].append(entry.as_list())

        if self.user_ccd is not None:
            content["userCCD"] = self.user_ccd

        return content

    @staticmethod
    def read(filename) -> InputFile:
        """
        Reads a file and returns its content as an instance of InputFile.

        Parameters
        ----------
        filename : str
            The path to the JSON file to be read.

        Returns
        -------
        InputFile
            An instance of InputFile that contains the deserialized content
            of the provided JSON file.
        """
        from .io import read_json
        return read_json(filename)

    def write(self, filename) -> None:
        """
        Writes the current object's data to a JSON file.

        Parameters
        ----------
        filename : str
            The path to the file where the object's data will be saved.
            The file must be writable.
        """
        from .io import write_json
        write_json(filename, self)
