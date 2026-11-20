from af3cli import Atom, Bond, InputFile


def test_default_seed_is_one() -> None:
    afinput = InputFile()
    assert afinput.seeds == {1}
    assert afinput.to_dict()["modelSeeds"] == [1]


def test_default_version_is_one() -> None:
    afinput = InputFile()
    assert afinput.version == 1
    assert afinput.to_dict()["version"] == 1


def test_merge_ignores_seeds_by_default() -> None:
    afinput = InputFile(seeds=[1])
    other = InputFile(seeds=[7, 8])
    afinput.merge(other)
    assert afinput.seeds == {1}


def test_merge_ignores_bonded_atoms_by_default() -> None:
    afinput = InputFile()
    other = InputFile()
    other.bonded_atoms.append(Bond(Atom("A", 1, "CA"), Atom("B", 2, "O")))
    afinput.merge(other)
    assert afinput.bonded_atoms == []
    assert "bondedAtomPairs" not in afinput.to_dict()
