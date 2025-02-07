from __future__ import annotations

from enum import StrEnum
from abc import ABCMeta
from typing import Generator
import re

from .mixin import DictMixin
from .exception import (AFSequenceError, AFTemplateError,
                        AFModificationError)
from .seqid import IDRecord


class SequenceType(StrEnum):
    """
    Represents the types of sequences.

    This enumeration defines constants for types of available sequences
    in the AlphaFold3 input file.
    """
    PROTEIN: str = "protein"
    RNA: str = "rna"
    DNA: str = "dna"


class TemplateType(StrEnum):
    """
    Represents both types of templates that can be used for
    protein sequences.

    Attributes
    ----------
    FILE : str
        Represents a file-based template with an absolute or relative
        path to the input file.
    STRING : str
        Represents a string-based template for inline use of the
        mmCIF file.
    """
    FILE: str = "mmcifPath"
    STRING: str = "mmcif"


class Template(DictMixin):
    """
    Manages structural templates for protein sequences. The `mmcif` attribute
    is either a string or a file path, depending on the template type.

    Attributes
    ----------
    template_type : TemplateType
        The type of template, indicating a string or file-based template.
    mmcif : str
        The mmCIF formatted string of the template or a file path.
    qidx : list of int
        List of query indices that correspond to specific positions in
        the query structure.
    tidx : list of int
        List of template indices that map to the positions in the template
        structure.
    """
    def __init__(
        self,
        template_type: TemplateType,
        mmcif: str,
        qidx: list[int],
        tidx: list[int]
    ):
        self.template_type: TemplateType = template_type
        self.mmcif: str = mmcif
        self.qidx: list[int] = qidx
        self.tidx: list[int] = tidx

    def to_dict(self):
        """
        Converts the attributes of the object into a dictionary representation
        to automatically generate the corresponding fields in the AlphaFold3 input.

        Returns
        -------
        dict
            A dictionary containing key-value pairs derived from the object's
            attributes. The keys can include 'mmcifPath', 'mmcif',
            'queryIndices', and 'templateIndices', depending on the values and
            conditions of the attributes.
        """
        return {
            self.template_type.value: self.mmcif,
            "queryIndices": self.qidx,
            "templateIndices": self.tidx
       }

    def __str__(self) -> str:
        return f"Template({self.template_type.value})"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}({self.template_type.name})>"


class MSA(DictMixin):
    """
    Manages paired and unpaired Multiple Sequence Alignment (MSA) data.
    These are only available for protein and RNA sequences and not for
    DNA sequences. A check will be performed in the `Sequence` class.

    Attributes
    ----------
    paired : str or None
        Paired MSA data or its file path.
    unpaired : str or None
        Unpaired MSA data or its file path.
    paired_is_path : bool
        Indicates whether `paired` represents a file path.
    unpaired_is_path : bool
        Indicates whether `unpaired` represents a file path.
    """
    def __init__(
        self,
        paired: str | None = None,
        unpaired: str | None = None,
        paired_is_path: bool = False,
        unpaired_is_path: bool = False,
    ):
        self.paired: str | None = paired
        self.unpaired: str | None = unpaired
        self.paired_is_path: bool = paired_is_path
        self.unpaired_is_path: bool = unpaired_is_path

    def to_dict(self) -> dict:
        """
        Converts the attributes of the object into a dictionary representation
        to automatically generate the corresponding fields in the AlphaFold3 input.

        Returns
        -------
        dict
            A dictionary containing key-value pairs derived from the object's
            attributes. The keys can include 'pairedMsaPath', 'pairedMsa',
            'unpairedMsaPath', or 'unpairedMsa', depending on the values and
            conditions of the attributes.
        """
        tmp_dict = {}
        if self.paired is not None:
            if self.paired_is_path:
                tmp_dict["pairedMsaPath"] = self.paired
            else:
                tmp_dict["pairedMsa"] = self.unpaired
        if self.unpaired is not None:
            if self.unpaired_is_path:
                tmp_dict["unpairedMsaPath"] = self.unpaired
            else:
                tmp_dict["unpairedMsa"] = self.unpaired
        return tmp_dict

    def __str__(self) -> str:
        display_paired = "paired" if self.paired is not None else ""
        display_unpaired = "unpaired" if self.unpaired is not None else ""
        return f"MSA({display_paired}, {display_unpaired})"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}>"


class Modification(DictMixin, metaclass=ABCMeta):
    """
    Represents a modification with specific CCD code and position.

    Attributes
    ----------
    mod_str : str
        The CCD code representing the type or name of the modification.
    mod_pos : int
        The position of the modification within its given context.
    """
    def __init__(self, mod_str: str, mod_pos: int):
        self.mod_str: str = mod_str
        self.mod_pos: int = mod_pos


class ResidueModification(Modification):
    """
    Represents a specific modification at a specific residue position of
    protein sequences.

    Attributes
    ----------
    mod_str : str
        The CCD code representing the type of the modification
        (e.g., phosphorylation, methylation).
    mod_pos : int
        An integer representing the position of the modification
        within the sequence.
    """
    def __init__(self, mod_str: str, mod_pos: int):
        super().__init__(mod_str, mod_pos)

    def to_dict(self):
        return {
            "ptmType": self.mod_str,
            "ptmPosition": self.mod_pos
        }

    def __str__(self) -> str:
        return f"ResidueModification({self.mod_str}, {self.mod_pos})"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}>"


class NucleotideModification(Modification):
    """
    Represents a specific modification at a specific residue position of
    nucleotide sequences.

    Attributes
    ----------
    mod_str : str
        Type of the nucleotide modification.
        :TODO: provide a list of possible modifications, a reference, or some common examples.
    mod_pos : int
        Position of the modification in the nucleotide sequence.
    """
    def __init__(self, mod_str: str, mod_pos: int):
        super().__init__(mod_str, mod_pos)

    def to_dict(self):
        return {
            "modificationType": self.mod_str,
            "basePosition": self.mod_pos
        }

    def __str__(self) -> str:
        return f"NucleotideModification({self.mod_str}, {self.mod_pos})"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}>"


class Sequence(IDRecord, DictMixin):
    """
    Represents a sequence with templates and modifications.

    The Sequence class is used to store any sequence data with associated
    modifications, templates, and multiple sequence alignment (MSA).
    It provides functionality for validation of sequence attributes,
    This class extends `_IDRecord` to automatically handle sequence IDs.

    Attributes
    ----------
    _seq_type : SequenceType
        The type of the sequence (e.g., Protein, DNA, RNA).
    _seq_str : str
        The string representation of the sequence.
    _msa : MSA or None
        The multiple sequence alignment (MSA) information, if available.
    _modifications : list of Modification
        Modifications associated with the sequence.
    _templates : list of Template
        Templates associated with the sequence. Supported only for protein sequences.
    _num : int
        The number of sequences associated with the sequence ID. This value
        will be overwritten if `seq_id` is specified.
    _seq_id : list[str] or None
        The sequence ID(s) associated with the sequence. These can be
        either specified as a list of strings or will be automatically
        assigned by `IDRegister`.
    """
    def __init__(
        self,
        seq_type: SequenceType,
        seq_str: str,
        seq_name: str = "",
        num: int = 1,
        seq_id: list[str] | None = None,
        modifications: list[Modification] | None = None,
        templates: list[Template] | None = None,
        msa: MSA | None = None,
    ):
        super().__init__(num, None)
        self._seq_str: str = seq_str
        self._seq_type: SequenceType = seq_type
        self.seq_name: str = seq_name

        self._msa: MSA | None = msa

        if modifications is None:
            modifications = []
        self._modifications: list[Modification] = modifications

        if seq_type != SequenceType.PROTEIN and templates is not None:
            raise AFTemplateError("Templates are only supported for proteins.")

        if templates is None:
            templates = []
        self._templates: list[Template] = templates

        if seq_id is not None:
            self._seq_id: list[str] = seq_id
            self._num: int = len(seq_id)

    @property
    def sequence(self) -> str:
        return self._seq_str

    @property
    def sequence_type(self) -> SequenceType:
        return self._seq_type

    @property
    def msa(self) -> MSA | None:
        return self._msa

    @property
    def modifications(self) -> list[Modification]:
        return self._modifications

    def _validate_modification_types(self):
        """
        Checks the validity of sequence modifications against the sequence type.

        Returns
        -------
        bool
            True if all modifications in `modifications` are valid for the given
            `seq_type`, or if `modifications` is empty. False otherwise.
        """
        if len(self._modifications) == 0:
            return True
        if self._seq_type == SequenceType.PROTEIN:
            return all(
                isinstance(mod, ResidueModification)
                for mod in self._modifications
            )
        else:
            return all(
                isinstance(mod, NucleotideModification)
                for mod in self._modifications
            )

    def to_dict(self) -> dict:
        """
         Convert the object to a dictionary representation.

         This method creates and returns a dictionary representation of the object,
         including its identifier, sequence string, modifications, templates, and,
         the multiple sequence alignment (MSA). It is used to generate the sequence
         entries in the AlphaFold3 input file.

         Returns
         -------
         dict
             A dictionary representation of the object.

         Raises
         ------
         AFSequenceError
             If the sequence is invalid for the sequence type.
         AFModificationError
             If the modifications are invalid for the sequence type.
         """
        if not is_valid_sequence(self._seq_type, self._seq_str):
            raise AFSequenceError(
                f"Invalid sequence for sequence type "
                f"{self._seq_type.name} ({self})."
            )

        if not self._validate_modification_types():
            raise AFModificationError(
                f"Invalid modification types for sequence {self}."
            )

        content = dict()
        content["id"] = self.get_id()
        content["sequence"] = self._seq_str
        if len(self._modifications):
            content["modifications"] = [m.to_dict() for m in self._modifications]
        if self._msa is not None:
            content |= self._msa.to_dict()
        return {self._seq_type.value: content}

    def __str__(self) -> str:
        display_template = "T" if len(self._templates) else ""
        display_mod = "M" if len(self._modifications) else ""
        display_msa = "MSA" if self._msa is not None else ""
        display_flags = ",".join(
            v for v in [display_template, display_mod, display_msa] if v
        )
        return f"{self._seq_type.name}({len(self._seq_str)})[{display_flags}]"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}({self._seq_type.name})>"


class ProteinSequence(Sequence):
    """
    Represents a protein sequence, inheriting properties and functionality
    from the `Sequence` base class.
    """
    def __init__(
        self,
        seq_str: str,
        seq_name: str | None = None,
        num: int = 1,
        seq_id: list[str] | None = None,
        modifications: list[ResidueModification] | None = None,
        templates: list[Template] | None = None,
        msa: MSA | None = None,
    ):
        super().__init__(
            SequenceType.PROTEIN,
            seq_str=seq_str,
            seq_name=seq_name,
            num=num,
            seq_id=seq_id,
            modifications=modifications,
            templates=templates,
            msa=msa,
        )

    @property
    def templates(self) -> list[Template]:
        return self._templates

    def to_dict(self) -> dict:
        content = super().to_dict().get(self._seq_type.value)
        if len(self._templates):
            content["templates"] = [t.to_dict() for t in self._templates]
        return {self._seq_type.value: content}


class DNASequence(Sequence):
    """
    Represents a DNA sequence, inheriting properties and functionality
    from the `Sequence` base class.
    """
    def __init__(
        self,
        seq_str: str,
        seq_name: str | None = None,
        num: int = 1,
        seq_id: list[str] | None = None,
        modifications: list[NucleotideModification] | None = None,
    ):
        super().__init__(
            SequenceType.DNA,
            seq_str=seq_str,
            seq_name=seq_name,
            num=num,
            seq_id=seq_id,
            modifications=modifications,
        )

    def reverse_complement(self) -> DNASequence:
        """
            Generates the reverse complement of the DNA sequence.

            Returns
            -------
            DNASequence
                A new DNASequence instance representing the reverse complement
                of this sequence with the same sequence name and number of
                entities.
        """
        cmap = {
            'A': 'T',
            'T': 'A',
            'C': 'G',
            'G': 'C'
        }
        complement = ''.join(cmap[base] for base in self._seq_str)
        complement = complement[::-1]
        return DNASequence(
            complement, num=self._num, seq_name=self.seq_name,
        )


class RNASequence(Sequence):
    """
    Represents a RNA sequence, inheriting properties and functionality
    from the `Sequence` base class.
    """
    def __init__(
        self,
        seq_str: str,
        seq_name: str | None = None,
        num: int = 1,
        seq_id: list[str] | None = None,
        modifications: list[NucleotideModification] | None = None,
        msa: MSA | None = None,
    ):
        super().__init__(
            SequenceType.RNA,
            seq_str=seq_str,
            seq_name=seq_name,
            num=num,
            seq_id=seq_id,
            modifications=modifications,
            msa=msa,
        )

    def reverse_complement(self) -> RNASequence:
        """
        Generates the reverse complement of the RNA sequence.

        Returns
        -------
        RNASequence
            A new RNASequence instance representing the reverse complement
            of this sequence with the same sequence name and number of
            entities.
        """
        cmap = {
            'A': 'U',
            'U': 'A',
            'C': 'G',
            'G': 'C'
        }
        complement = ''.join(cmap[base] for base in self._seq_str)
        complement = complement[::-1]
        return RNASequence(
            complement, num=self._num, seq_name=self.seq_name,
        )


def read_fasta(filename: str) -> Generator[tuple[str, str], None, None]:
    """
    Reads a FASTA file and yields sequences with their identifiers as tuples.

    This function utilizes Biopython's SeqIO to parse a given FASTA file and
    yield tuples containing sequence identifiers and their corresponding
    sequences. The function ensures that Biopython is installed and raises
    an appropriate error if it is not available.

    Parameters
    ----------
    filename : str
        Path to the FASTA file to be read.

    Yields
    ------
    tuple of (str, str)
        A tuple where the first element is the sequence identifier, and the
        second element is the sequence as a string.

    Raises
    ------
    ImportError
        If Biopython is not installed and cannot be imported.
    """
    try:
        from Bio import SeqIO
        for entry in SeqIO.parse(filename, "fasta"):
            yield entry.id, str(entry.seq).upper()
    except ImportError as e:
        raise ImportError("Please install Biopython to read FASTA files") from e


def is_valid_sequence(seq_type: SequenceType, seq_str: str) -> bool:
    """
    Determines if a given sequence string corresponds to the specified sequence type.

    Parameters
    ----------
    seq_type : SequenceType
        Type of the sequence to validate.
    seq_str : str
        The sequence string to validate against the specified sequence type.

    Returns
    -------
    bool
        True if all characters in `seq_str` belong to the valid character set for
        the specified `seq_type`; False otherwise.
    """
    SEQ_CHAR_SETS = {
        SequenceType.PROTEIN: set("ACDEFGHIKLMNPQRSTVWY"),
        SequenceType.DNA: set("ACGT"),
        SequenceType.RNA: set("ACGU")
    }
    return all(char in SEQ_CHAR_SETS[seq_type] for char in seq_str)


def identify_sequence_type(seq_str: str) -> SequenceType | None:
    """
    Identifies the type of a given biological sequence based on its composition.
    The function examines if the sequence can be classified as DNA, RNA, or
    protein, and returns the corresponding sequence type.

    Parameters
    ----------
    seq_str : str
        The biological sequence to be analyzed.

    Returns
    -------
    SequenceType or None
        The type of the sequence, identified as one of the following:
        - SequenceType.DNA: If the sequence is identified as DNA.
        - SequenceType.RNA: If the sequence is identified as RNA.
        - SequenceType.PROTEIN: If the sequence is identified as protein.
        Returns None if the sequence is ambiguous (e.g., qualifies as both DNA and RNA)
        or does not fit any of the known sequence types.
    """
    is_protein = is_valid_sequence(SequenceType.PROTEIN, seq_str)
    is_dna = is_valid_sequence(SequenceType.DNA, seq_str)
    is_rna = is_valid_sequence(SequenceType.RNA, seq_str)

    if is_dna and is_rna:
        return None
    elif is_dna:
        return SequenceType.DNA
    elif is_rna:
        return SequenceType.RNA
    elif is_protein:
        return SequenceType.PROTEIN
    return None


def sanitize_sequence_name(name: str) -> str:
    """
    Sanitizes a sequence name by replacing unwanted characters and
    stripping whitespace.

    Parameters
    ----------
    name : str
        The original sequence name to be sanitized.

    Returns
    -------
    str
        The sanitized sequence name.
    """
    return re.sub(r'[ |:|(|)|]', '_', name).strip()

def fasta2seq(filename: str) -> Generator[Sequence | None, None, None]:
    """
    Converts a FASTA file into a sequence generator.

    Parameters
    ----------
    filename : str
        The path to the FASTA file to be read.

    Yields
    ------
    Sequence or None
        A `Sequence` object if the sequence type can be identified; otherwise, `None`.
    """
    for entry_name, entry_seq in read_fasta(filename):
        if entry_seq is None:
            yield None
            continue

        seq_type = identify_sequence_type(entry_seq)
        if seq_type is None:
            yield None
            continue        

        match seq_type:
            case SequenceType.PROTEIN:
                yield ProteinSequence(seq_str=entry_seq, seq_name=sanitize_sequence_name(entry_name))
            case SequenceType.DNA:
                yield DNASequence(seq_str=entry_seq, seq_name=sanitize_sequence_name(entry_name))
            case SequenceType.RNA:
                yield RNASequence(seq_str=entry_seq, seq_name=sanitize_sequence_name(entry_name))
        

def read_first_seq_fasta(filename: str, num: int = 1) -> Sequence:
    """
    Converts a FASTA file into a sequence generator.

    Parameters
    ----------
    filename : str
        The path to the FASTA file to be read.

    Yields
    ------
    Sequence or None
        A `Sequence` object if the sequence type can be identified; otherwise, `None`.
    """
    
    seq_tuple = next(read_fasta(filename))

    seq_type = identify_sequence_type(seq_tuple[1])
    if seq_type is None:
        return None

    match seq_type:
        case SequenceType.PROTEIN:
            return ProteinSequence(seq_str=seq_tuple[1], seq_name=sanitize_sequence_name(seq_tuple[0]), num=num)
        case SequenceType.DNA:
            return DNASequence(seq_str=seq_tuple[1], seq_name=sanitize_sequence_name(seq_tuple[0]), num=num)
        case SequenceType.RNA:
            return RNASequence(seq_str=seq_tuple[1], seq_name=sanitize_sequence_name(seq_tuple[0]), num=num)