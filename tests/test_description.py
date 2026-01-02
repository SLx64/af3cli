import pytest

from af3cli import InputFile, ProteinSequence, SMILigand


@pytest.fixture
def descriptions() -> tuple[str, str]:
    return "Protein description", "Ligand description"


@pytest.mark.parametrize(
    "version, expect_description",
    [
        (None, False),   # default case (version implicitly set to 1)
        (2, False),      # version == 2
        (4, True),       # version == 4
    ],
)
def test_description_handling(version: int | None, expect_description: bool, descriptions: tuple[str, str]) -> None:
    seq_desc, lig_desc = descriptions

    afinput = InputFile() if version is None else InputFile(version=version)
    afinput.sequences.append(ProteinSequence("ACDEFGHIK", description=seq_desc))
    afinput.ligands.append(SMILigand("CCO", description=lig_desc))

    content = afinput.to_dict()

    # Extract the inner entries from sequences
    inners = [entry for sequence in content.get("sequences", []) for entry in sequence.values()]

    if expect_description:
        assert seq_desc in [inner.get("description") for inner in inners]
        assert lig_desc in [inner.get("description") for inner in inners]
    else:
        assert all("description" not in inner for inner in inners)
