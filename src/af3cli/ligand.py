from enum import StrEnum

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

    The Ligand class is used to represent a specific ligand with its type,
    string representation, sequence ID, and count.

    Attributes
    ----------
    ligand_str : list of str or str
        The string representation(s) of the ligand.
    ligand_type : LigandType
        The type of the ligand entry.
    num : int or None
        The number of ligand sequences, default is 1. If `seq_id` is provided and
        `num` is not, then `num` will be inferred from the length of `seq_id`.
    _seq_id : list[str] or None
        The sequence ID(s) associated with the sequence. These can be
        either specified as a list of strings or will be automatically
        assigned by `IDRegister`.
    """
    def __init__(
        self,
        ligand_type: LigandType,
        ligand_str: list[str] | str,
        num: int | None = None,
        seq_id: list[str] | None = None
    ):
        super().__init__(None)
        self.ligand_str: list[str] | str  = ligand_str
        self.ligand_type: LigandType = ligand_type

        # can be overwritten if seq_id is specified
        if num is None:
            self.num: int = 1
        else:
            self.num: int = num

        if seq_id is not None:
            self._seq_id: list[str] = seq_id
            if num is None:
                self.num: int = len(seq_id)
            elif len(seq_id) != num:
                raise ValueError((f"Sequence ID length ({len(seq_id)}) does "
                                  f"not match sequence number ({num})."))

    def to_dict(self):
        """
        Converts the object's data into a dictionary format to automatically
        convert it to the AlphaFold3 input file format.

        Returns
        -------
        dict
            A dictionary containing the object's ID and ligand data.
        """
        content = dict()
        content["id"] = self.get_id()
        content[self.ligand_type.value] = self.ligand_str
        return {"ligand": content}
