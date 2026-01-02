import pytest

import io

from af3cli.sequence import (SequenceType, Sequence,
                             ProteinSequence,
                             DNASequence, RNASequence)
from af3cli.sequence import TemplateType, Template
from af3cli.sequence import (Modification, ResidueModification,
                             NucleotideModification)
from af3cli.sequence import MSA
from af3cli.sequence import is_valid_sequence, identify_sequence_type
from af3cli.sequence import read_fasta, fasta2seq


@pytest.fixture(scope="module")
def single_fasta_file():
    content = (">sp|Q9C0K0|BC11B_HUMAN B-cell\n"
               "MSRRKQGNPQHLSQRELITPEADH")

    return io.StringIO(content)


@pytest.fixture(scope="module")
def multi_protein_fasta_file():
    content = (">sp|Q9H165|BC11A_HUMAN\n"
               "MSRRKQGKPQHLSKREFSPEPLEA\n"
               ">sp|Q9C0K0|BC11B_HUMAN\n"
               "MSRRKQGNPQHLSQRELITPEADH")
    return io.StringIO(content)


@pytest.fixture(scope="module")
def multi_all_fasta_file():
    content = (">protein\n"
               "MSRRKQGKPQHLSKREFSPEPLEA\n"
               ">dna\n"
               "AATTTTCCCGGGGGTTTTTTAAAA\n"
               ">rna\n"
               "AACCCGGGGUUUGGCCGGGAAUUU")
    return io.StringIO(content)


@pytest.mark.parametrize("seq_type,seq_type_value", [
    (SequenceType.PROTEIN, "protein"),
    (SequenceType.RNA, "rna"),
    (SequenceType.DNA, "dna")
])
def test_residue_type(seq_type: SequenceType, seq_type_value: str) -> None:
    assert SequenceType(seq_type).value == seq_type_value


@pytest.mark.parametrize("template_type,template_type_value",[
    (TemplateType.STRING, "mmcif"),
    (TemplateType.FILE, "mmcifPath")
])
def test_template_type(
        template_type: TemplateType,
        template_type_value: str
) -> None:
    assert TemplateType(template_type).value == template_type_value


@pytest.mark.parametrize("template_type,mmcif,qidx,tidx", [
    (TemplateType.STRING, "data_ ...", [1, 2, 3, 4], [4, 6, 7, 8]),
    (TemplateType.FILE, "/path_to_file", [1, 2, 3, 4], [4, 6, 7, 8]),
])
def test_template_init(
        template_type: TemplateType,
        mmcif: str,
        qidx: list[int],
        tidx: list[int]
) -> None:
    template = Template(template_type, mmcif, qidx, tidx)
    assert template.template_type == template_type
    assert template.mmcif == mmcif
    assert template.qidx == qidx
    assert template.tidx == tidx


@pytest.mark.parametrize("template_type,mmcif,qidx,tidx", [
    (TemplateType.STRING, "data_ ...", [1, 2, 3, 4], [4, 6, 7, 8]),
    (TemplateType.FILE, "/path_to_file", [1, 2, 3, 4], [4, 6, 7, 8]),
])
def test_template_to_dict(
        template_type: TemplateType,
        mmcif: str,
        qidx: list[int],
        tidx: list[int]
) -> None:
    template = Template(template_type, mmcif, qidx, tidx)
    tdict = template.to_dict()
    assert isinstance(tdict, dict)
    assert template_type.value in tdict.keys()
    assert "queryIndices" in tdict.keys()
    assert "templateIndices" in tdict.keys()
    assert tdict[template_type.value] == mmcif
    assert tdict["queryIndices"] == qidx
    assert tdict["templateIndices"] == tidx


@pytest.mark.parametrize("mtype,mpos", [
    ("HY3", 1), ("P1L", 5)
])
def test_residue_mod(mtype: str, mpos: int) -> None:
    modification = ResidueModification(mtype, mpos)
    assert modification.mod_str == mtype
    assert modification.mod_pos == mpos


@pytest.mark.parametrize("mtype,mpos", [
    ("HY3", 1), ("P1L", 5)
])
def test_residue_mod_dict(mtype: str, mpos: int) -> None:
    modification = ResidueModification(mtype, mpos)
    mdict = modification.to_dict()
    assert isinstance(mdict, dict)
    assert "ptmType" in mdict.keys()
    assert "ptmPosition" in mdict.keys()
    assert mdict["ptmType"] == mtype
    assert mdict["ptmPosition"] == mpos


@pytest.mark.parametrize("mtype,mpos", [
    ("6OG", 1), ("6MA", 2), ("2MG", 5), ("5MC", 10)
])
def test_nucleotide_mod(mtype: str, mpos: int) -> None:
    modification = NucleotideModification(mtype, mpos)
    assert modification.mod_str == mtype
    assert modification.mod_pos == mpos


@pytest.mark.parametrize("mtype,mpos", [
    ("6OG", 1), ("6MA", 2), ("2MG", 5), ("5MC", 10)
])
def test_nucleotide_mod_dict(mtype: str, mpos: int) -> None:
    modification = NucleotideModification(mtype, mpos)
    mdict = modification.to_dict()
    assert isinstance(mdict, dict)
    assert "modificationType" in mdict.keys()
    assert "basePosition" in mdict.keys()
    assert mdict["modificationType"] == mtype
    assert mdict["basePosition"] == mpos


@pytest.mark.parametrize("seq_type,seq_str,seq_id,seq_mods",[
    (SequenceType.PROTEIN, "MVKVGVNGF", ["A", "B"],
     [ResidueModification("HY3", 1),
      ResidueModification("P1L", 5)]),
    (SequenceType.RNA, "AUGUGUAU", ["A", "B"], []),
    (SequenceType.DNA, "GACCTCT", None,
     [NucleotideModification("6OG", 1),
      NucleotideModification("6MA", 2)])
])
def test_sequence_init_mod(
        seq_type: SequenceType,
        seq_str: str, seq_id: list[str] | str | None,
        seq_mods: list[Modification]
) -> None:
    seq = Sequence(seq_type, seq_str, seq_id=seq_id, modifications=seq_mods)
    assert seq._seq_type == seq_type
    assert seq._seq_str == seq_str
    assert seq._modifications == seq_mods
    if isinstance(seq_id, list) or seq_id is None:
        assert seq.get_id() == seq_id
    else:
        assert seq.get_id() == [seq_id]


@pytest.mark.parametrize("seq_type,seq_str,seq_id,templates",[
    (SequenceType.PROTEIN, "MVKVGVNGF", ["A", "B"],
     [Template(TemplateType.STRING, "data", [1, 2, 3, 4], [4, 6, 7, 8]),
      Template(TemplateType.FILE, "/path", [1, 2, 3, 4], [4, 6, 7, 8])]
    ),
    (SequenceType.PROTEIN, "MVKVGVNGF", ["A", "B"],
     [Template(TemplateType.STRING, "data", [1, 2], [4, 6])]
    ),
])
def test_sequence_init_template(
        seq_type: SequenceType,
        seq_str: str,
        seq_id: list[str] | str | None,
        templates: list[Template] | None
) -> None:
    seq = Sequence(seq_type, seq_str, seq_id=seq_id, templates=templates)
    assert seq.sequence_type == seq_type
    assert seq.sequence == seq_str
    if isinstance(seq_id, list) or seq_id is None:
        assert seq.get_id() == seq_id
    else:
        assert seq.get_id() == [seq_id]
    assert len(seq._templates) == len(templates)


@pytest.mark.parametrize("paired,unpaired,pispath,unpispath", [
    (None, None, False, False),
    ("test", None, False, False),
    (None, "test", False, False),
    ("test", "test", True, False),
    ("test", "test", False, True),
    ("test", "test", True, True),
])
def test_msa_init(
        paired: str | None,
        unpaired: str | None,
        pispath: bool,
        unpispath: bool
) -> None:
    msa = MSA(paired, unpaired, pispath, unpispath)
    assert msa.paired == paired
    assert msa.unpaired == unpaired
    assert msa.paired_is_path == pispath
    assert msa.unpaired_is_path == unpispath


@pytest.mark.parametrize("paired,unpaired,pispath,unpispath", [
    (None, None, False, False),
    ("test", None, False, False),
    (None, "test", False, False),
    ("test", "test", True, False),
    ("test", "test", False, True),
    ("test", "test", True, True),
])
def test_msa_to_dict(
        paired: str | None,
        unpaired: str | None,
        pispath: bool,
        unpispath: bool
) -> None:
    msa = MSA(paired, unpaired, pispath, unpispath)
    tmp_dict = msa.to_dict()
    if paired is not None:
        if pispath:
            assert "pairedMsaPath" in tmp_dict.keys()
            assert "pairedMsa" not in tmp_dict.keys()
        else:
            assert "pairedMsaPath" not in tmp_dict.keys()
            assert "pairedMsa" in tmp_dict.keys()

    if unpaired is not None:
        if unpispath:
            assert "unpairedMsaPath" in tmp_dict.keys()
            assert "unpairedMsa" not in tmp_dict.keys()
        else:
            assert "unpairedMsaPath" not in tmp_dict.keys()
            assert "unpairedMsa" in tmp_dict.keys()


@pytest.mark.parametrize("seq_type,seq_str,seq_id,msa",[
    (SequenceType.PROTEIN, "MVKVGVNGF", "A", None),
    (SequenceType.PROTEIN, "MVKVGVNGF", ["A", "B"],
     MSA(paired="test")
    ),
    (SequenceType.PROTEIN, "MVKVGVNGF", "A",
     MSA(paired="test", paired_is_path=True)
    ),
    (SequenceType.PROTEIN, "MVKVGVNGF", ["A", "B"],
     MSA(paired="test", unpaired="test")
    ),
    (SequenceType.RNA, "AUGUGUAU", "C",
     MSA(unpaired="test", unpaired_is_path=True)
    ),
    (SequenceType.PROTEIN, "MVKVGVNGF", ["A"],
     MSA(paired="test", unpaired="test", unpaired_is_path=True)
    ),
    (SequenceType.PROTEIN, "MVKVGVNGF", ["A"],
     MSA(paired="test", unpaired="test", paired_is_path=True, unpaired_is_path=True)
    ),
])
def test_sequence_init_msa(
        seq_type: SequenceType,
        seq_str: str,
        seq_id: list[str] | str | None,
        msa: MSA | None
) -> None:
    seq = Sequence(seq_type, seq_str, seq_id=seq_id, msa=msa)
    assert seq.sequence_type == seq_type
    assert seq.sequence == seq_str
    if isinstance(seq_id, list) or seq_id is None:
        assert seq.get_id() == seq_id
    else:
        assert seq.get_id() == [seq_id]
    tmp_dict = seq.to_dict()[seq_type.value]
    if msa is None:
        assert seq._msa is None
        assert "pairedMsa" not in tmp_dict.keys()
        assert "unpairedMsa" not in tmp_dict.keys()
        return
    if msa.paired_is_path:
        assert "pairedMsaPath" in tmp_dict.keys()
        assert "pairedMsa" not in tmp_dict.keys()
    if msa.unpaired_is_path:
        assert "unpairedMsaPath" in tmp_dict.keys()
        assert "unpairedMsa" not in tmp_dict.keys()

@pytest.mark.parametrize("seq_type,seq_str,seq_id", [
    (SequenceType.PROTEIN, "MVKVGVNGF", None),
    (SequenceType.PROTEIN, "AAQAA", ["A", "B"]),
    (SequenceType.RNA, "AUGUGUAU", "A"),
    (SequenceType.DNA, "GACCTCT", ["C"])
])
def test_sequence_remove_id(
        seq_type: SequenceType,
        seq_str: str,
        seq_id: list[str] | str | None,
):
    seq = Sequence(seq_type, seq_str, seq_id=seq_id)
    seq.remove_id()
    assert seq.get_id() is None

@pytest.mark.parametrize("seq_type,seq_str,result", [
    (SequenceType.PROTEIN, "MVKVGVNGF", True),
    (SequenceType.RNA, "AUGUGUAU", True),
    (SequenceType.DNA, "GACCTCT", True),
    (SequenceType.PROTEIN, "AAQAAU", False),
    (SequenceType.RNA, "ABCDEFG", False),
    (SequenceType.DNA, "1234567", False),
])
def test_is_seq_type(seq_type: SequenceType, seq_str: str, result: bool):
    assert is_valid_sequence(seq_type=seq_type, seq_str=seq_str) == result


@pytest.mark.parametrize("seq_str,result", [
    ("MVKVGVNGF", SequenceType.PROTEIN),
    ("AUGUGUAU", SequenceType.RNA),
    ("GACCTCT", SequenceType.DNA),
    ("GACCCAAGG", None),
    ("1234567", None),
    ("ABCDEFG", None),
])
def test_identify_seq_type(seq_str: str, result: SequenceType | None):
    assert identify_sequence_type(seq_str=seq_str) == result


def test_read_fasta_single(single_fasta_file):
    for entry in read_fasta(single_fasta_file):
        assert entry[0] is not None
        assert entry[1] is not None
        assert isinstance(entry[0], str)
        assert isinstance(entry[1], str)
        assert entry[1] == "MSRRKQGNPQHLSQRELITPEADH"


def test_fasta2seq_single(single_fasta_file):
    for entry in fasta2seq(single_fasta_file):
        assert entry is not None
        assert isinstance(entry, Sequence)
        assert entry.sequence_type == SequenceType.PROTEIN
        assert entry.sequence == "MSRRKQGNPQHLSQRELITPEADH"


def test_fasta2seq_multi(multi_protein_fasta_file):
    for entry in fasta2seq(multi_protein_fasta_file):
        assert entry is not None
        assert isinstance(entry, Sequence)
        assert entry.sequence_type == SequenceType.PROTEIN


def test_fasta2seq_multi_all(multi_all_fasta_file):
    for entry in fasta2seq(multi_all_fasta_file):
        assert entry is not None
        assert isinstance(entry, Sequence)
        assert entry.sequence_type is not None


@pytest.mark.parametrize("cls,seq_str,seq_type", [
    (ProteinSequence, "MVKVGVNGF", SequenceType.PROTEIN),
    (DNASequence, "AUGUGUAU", SequenceType.DNA),
    (RNASequence, "GACCTCT", SequenceType.RNA)
])
def test_sequence_shortcuts(cls, seq_str: str, seq_type: SequenceType):
    seq = cls(seq_str)
    assert seq.sequence_type == seq_type
    assert seq.sequence == seq_str
    assert seq.get_id() is None
    assert len(seq.modifications) == 0
    assert seq.msa is None
    assert seq.num == 1
    if isinstance(seq, ProteinSequence):
        assert len(seq.templates) == 0


def test_dna_complement():
    seq = DNASequence("GCGAATTCG")
    complement = seq.reverse_complement()
    assert complement.sequence == "CGAATTCGC"


@pytest.mark.parametrize("cls,seq_str,desc", [
    (ProteinSequence, "MVKVGVNGF", "kinase protein"),
    (DNASequence, "GACCTCT", "plasmid insert"),
    (RNASequence, "AUGUGUAU", None),
])
def test_sequence_description_to_dict(cls, seq_str: str, desc: str | None):
    seq = cls(seq_str=seq_str, description=desc)
    ddict = seq.to_dict()
    key = seq.sequence_type.value
    assert key in ddict
    inner = ddict[key]
    assert inner["sequence"] == seq_str
    if desc:
        assert inner["description"] == desc



@pytest.mark.parametrize("cls,seq_str,seq_id,num", [
    (ProteinSequence, "MVKVGVNGF", "A", 1),
    (DNASequence, "AUGUGUAU", ["A", "B"], 1),
    (RNASequence, "GACCTCT", None, 4)
])
def test_num_sequences(
        cls: type[ProteinSequence | DNASequence | RNASequence],
        seq_str: str, seq_id: str | list[str], num: int):
    seq = cls(seq_str=seq_str, seq_id=seq_id, num=num)
    if seq_id is None:
        assert seq.num == num
    else:
        num_ids = len(seq.get_id())
        if num_ids > num:
            assert seq.num == num_ids
        else:
            assert seq.num == num
