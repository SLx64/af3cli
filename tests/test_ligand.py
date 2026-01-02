import pytest

from af3cli.ligand import Ligand, LigandType, CCDLigand, SMILigand


@pytest.mark.parametrize("lig_type,lig_type_value",[
    (LigandType.CCD, "ccdCodes"),
    (LigandType.SMILES, "smiles")
])
def test_ligand_type(lig_type: LigandType, lig_type_value: str) -> None:
    assert LigandType(lig_type).value == lig_type_value


@pytest.mark.parametrize("lig_type,lig_str,num,seq_id,actual_num",[
    (LigandType.CCD, ["NAC"], 1, None, 1),
    (LigandType.CCD, ["ATP"], 2, ["A", "B"], 2),
    (LigandType.SMILES, "CCC", 1, None, 1),
    (LigandType.SMILES, "CCC", 1, None, 1),
    (LigandType.SMILES, "CCC", 2, None, 2),
    (LigandType.SMILES, "CCC", 2, ["A", "B"], 2),
    (LigandType.SMILES, "CCC", 1, ["A", "B"], 2)
])
def test_ligand_init(
        lig_type: LigandType,
        lig_str: list[str] | str,
        num: int,
        seq_id: list[str] | str | None,
        actual_num: int
) -> None:
    ligand = Ligand(lig_type, lig_str, num, seq_id)
    assert ligand.ligand_type == lig_type
    assert ligand.ligand_value == lig_str
    assert ligand.num == actual_num
    assert ligand.get_id() == seq_id


@pytest.mark.parametrize("lig_type,lig_str,seq_id,desc",[
    (LigandType.CCD, ["ATP"], ["A", "B"], "energy molecule"),
    (LigandType.SMILES, "CCC", ["A", "B"], None),
    (LigandType.SMILES, "CCC", ["A", "B"], "short alkane")
])
def test_ligand_to_dict(
        lig_type: LigandType,
        lig_str: str,
        seq_id: list[str] | str | None,
        desc: str | None
) -> None:
    ligand = Ligand(lig_type, lig_str, description=desc, seq_id=seq_id)
    lig_dict = ligand.to_dict()
    key, values = next(iter(lig_dict.items()))
    assert key == "ligand"
    assert values[lig_type.value] == lig_str
    assert values["id"] == seq_id
    if desc:
        assert values["description"] == desc


@pytest.mark.parametrize("cls,lig_type,lig_str,num,seq_id,actual_num",[
    (CCDLigand, LigandType.CCD, ["NAC"], 1, None, 1),
    (CCDLigand, LigandType.CCD, ["CMPD1", "CMPD2"], 1, None, 1),
    (CCDLigand, LigandType.CCD, "NAC", 1, None, 1),
    (CCDLigand, LigandType.CCD, ["ATP"], 2, ["A", "B"], 2),
    (SMILigand, LigandType.SMILES, "CCC", 1, None, 1),
    (SMILigand, LigandType.SMILES, "CCC", 1, None, 1),
    (SMILigand, LigandType.SMILES, "CCC", 2, None, 2),
    (SMILigand, LigandType.SMILES, "CCC", 2, ["A", "B"], 2),
    (SMILigand, LigandType.SMILES, ["CCC"], 1, ["A", "B"], 2)
])
def test_ligand_child_classes(
        cls,
        lig_type: LigandType,
        lig_str: str,
        num: int,
        seq_id: list[str] | str | None,
        actual_num: int
) -> None:
    ligand = cls(lig_str, num, seq_id)
    assert ligand.ligand_type == lig_type
    assert ligand.ligand_value == lig_str
    assert ligand.num == actual_num
    assert ligand.get_id() == seq_id


def test_invalid_ligand() -> None:
    with pytest.raises(ValueError):
        ligand = SMILigand(["CCC", "CCO"])
        ligand.to_dict()
