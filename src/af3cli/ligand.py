from enum import StrEnum
from typing import Generator

from .mixin import DictMixin
from .seqid import IDRecord


class LigandType(StrEnum):
    """
    Enumeration of possible ligand types.

    Attributes
    ----------
    SMILES : str
        Ligand representation using SMILES string notation.
    CCD: str
        Ligand representation using CCD codes.
    """
    SMILES: str = "smiles"
    CCD: str = "ccdCodes"


class Ligand(IDRecord, DictMixin):
    """
    Represents a ligand with associated type, sequence ID, and other attributes.

    The Ligand class is used as a base class to represent a specific
    ligand with its type, string representation, sequence ID, and count.

    Attributes
    ----------
    _ligand_value : list of str or str
        The string representation(s) of the ligand.
    _ligand_value : LigandType
        The type of the ligand entry.
    description : str or None
        Optional free-text description of the ligand entry.
    _seq_id : list[str] or None
        The sequence ID(s) associated with the sequence. These can be
        either specified as a list of strings or will be automatically
        assigned by `IDRegister`.
    """
    def __init__(
        self,
        ligand_type: LigandType,
        ligand_value: list[str] | str,
        num: int = 1,
        seq_id: list[str] | None = None,
        description: str | None = None,
    ):
        super().__init__(num, None)
        self._ligand_value: list[str] | str  = ligand_value
        self._ligand_type: LigandType = ligand_type
        self.description: str | None = description

        # can be overwritten if length of seq_id is larger
        self.num = num
        self.set_id(seq_id)

    @property
    def ligand_type(self) -> LigandType:
        return self._ligand_type

    @property
    def ligand_value(self) -> list[str] | str:
        return self._ligand_value

    def to_dict(self):
        """
        Converts the object's data into a dictionary format to automatically
        convert it to the AlphaFold3 input file format.

        Returns
        -------
        dict
            A dictionary containing the object's ID and ligand data.
        """
        match self._ligand_type:
            case LigandType.CCD:
                if isinstance(self._ligand_value, str):
                    # otherwise the CCD name string will be treated as list of chars
                    self._ligand_value = [self._ligand_value]
            case LigandType.SMILES:
                if isinstance(self._ligand_value, list):
                    if len(self._ligand_value) != 1:
                        raise ValueError("SMILES value must be a single string")
                    self._ligand_value, _ = self._ligand_value
            case _:
                raise ValueError(f"Invalid ligand type: {self._ligand_type}")

        content = dict()
        content["id"] = self.get_full_id_list()
        content[self._ligand_type.value] = self._ligand_value
        if self.description:
            content["description"] = self.description
        return {"ligand": content}

    def __str__(self) -> str:
        return f"Ligand({self._ligand_type.name})"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}>"


class CCDLigand(Ligand):
    """
    Represents a CCD (Chemical Component Dictionary) ligand.
    """
    def __init__(
        self,
        ligand_value: list[str],
        num: int = 1,
        seq_id: list[str] | None = None,
        description: str | None = None,
    ):
        super().__init__(LigandType.CCD, ligand_value, num, seq_id, description=description)


class SMILigand(Ligand):
    """
    Represents a ligand that uses SMILES notation to define the chemical structure.
    """
    def __init__(
        self,
        ligand_value: str,
        num: int = 1,
        seq_id: list[str] | None = None,
        description: str | None = None,
    ):
        super().__init__(LigandType.SMILES, ligand_value, num, seq_id, description=description)


def sdf2smiles(filename: str) -> Generator[str | None, None, None]:
    """
    Reads a Structure Data File (SDF) and converts the molecules into SMILES format.

    This function uses RDKit to process the molecules in an SDF file and converts
    them into the SMILES string representation. If any molecule cannot be read
    from the file, it will be skipped with a warning message. The total number of
    successfully converted molecules will be logged. If RDKit is not installed,
    the function will terminate the program with an informative error message.

    Parameters
    ----------
    filename : str
        The path to the SDF file that needs to be read.

    Returns
    -------
    Generator of str or None
        A generator that yields SMILES strings representing the molecules
        contained in the specified SDF file.

    Raises
    ------
    ImportError
        If RDKit is not installed on the system.
    """
    try:
        from rdkit import Chem
        supplier = Chem.SDMolSupplier(filename)
        for mol in supplier:
            if mol is None:
                yield None
            yield Chem.MolToSmiles(mol)
    except ImportError as e:
        raise ImportError("Please install RDKit to read SDF files") from e
