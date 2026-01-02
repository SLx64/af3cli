import json

from .input import InputFile
from .bond import Atom, Bond
from .exception import AFMissingFieldError, AFTemplateError, AFMSAError
from .builder import InputBuilder
from .ligand import Ligand, SMILigand, CCDLigand
from .sequence import (Sequence,
                       ResidueModification,
                       NucleotideModification,
                       SequenceType, TemplateType,
                       Template, MSA, Modification)


def write_json(filename: str, data: InputFile) -> None:
    """
    Writes the contents of an InputFile object to the specified JSON file. The
    data is serialized into a JSON structure with readable indentation and then
    saved to a file.

    Parameters
    ----------
    filename : str
        The name of the file where the JSON data will be written.
    data : InputFile
        The InputFile object containing the data to be serialized and written
        to the JSON file.
    """
    with open(filename, 'w') as json_file:
        json.dump(data.to_dict(), json_file, indent=4)


def _read(filename: str) -> dict:
    """
    Reads a JSON file and returns its content as a dictionary.

    Parameters
    ----------
    filename : str
        The path to the JSON file that should be read.

    Returns
    -------
    dict
        A dictionary containing the parsed data from the JSON file.
    """
    with open(filename, "r") as json_file:
        data = json.load(json_file)
        return data


def _check_data(data: dict) -> None:
    """
    Checks the integrity of input data dictionary by verifying the presence of
    mandatory keys in AlphaFold3 input files.

    Parameters
    ----------
    data : dict
        The input data dictionary to be validated. It must contain the keys
        "name", "version", "dialect", "modelSeeds", and "sequences".

    Raises
    ------
    AFMissingFieldError
        If mandatory keys are not present in the input dictionary.
    """
    for key in ["name", "version", "dialect",
                "modelSeeds", "sequences"]:
        if key not in data:
            raise AFMissingFieldError(
                f"Missing field '{key}' in input file."
            )
    if len(data["sequences"]) == 0:
        raise AFMissingFieldError(
            f"Missing sequences in input file."
        )


def _get_id(seq_content: dict) -> list[str]:
    """
    Retrieve the identifier(s) from the given sequence content.

    This method retrieves the "id" key from the sequence content and ensures
    that it is returned as a list. If the "id" is already a list, it is returned
    unchanged. If the "id" is a single string, it is converted into a single-item
    list and returned.

    Parameters
    ----------
    seq_content : dict
        A dictionary containing sequence data with the "id" key.

    Returns
    -------
    list of str
        A list of sequence identifiers. If the "id" is a string, it is wrapped
        in a list. If already a list, it is returned directly.
    """
    seq_id = seq_content["id"]
    if isinstance(seq_id, str):
        seq_id = [seq_id]
    return seq_id


def _parse_ligand(seq_content: dict) -> Ligand:
    """
    Parses a ligand from the provided sequence content.

    This static method extracts data from a dictionary to create a `Ligand`
    object. It determines the ligand type based on the presence of specific
    keys ("ccdCodes" or "smiles") in the input dictionary. If both keys are
    missing, an exception is raised to indicate incomplete or invalid data.

    Parameters
    ----------
    seq_content : dict
        A dictionary containing the ligand data.

    Returns
    -------
    Ligand
        An instance of the `Ligand` class, initialized with specific ligand
        type, ligand string, and IDs.

    Raises
    ------
    AFMissingFieldError
        If neither "ccdCodes" nor "smiles" are present in the input dictionary.
    """
    seq_id = _get_id(seq_content)

    description = seq_content.get("description")

    if "ccdCodes" in seq_content.keys():
        return CCDLigand(
            ligand_value=seq_content["ccdCodes"],
            seq_id=seq_id,
            description=description
        )
    elif "smiles" in seq_content.keys():
        return SMILigand(
            ligand_value=seq_content["smiles"],
            seq_id=seq_id,
            description=description
        )
    else:
        raise AFMissingFieldError(
            f"Missing 'ccdCodes' or 'smiles' in sequence {seq_id}."
        )


def _parse_modifcations(
    seq_type: SequenceType,
    seq_content: dict
) -> list[Modification]:
    """
    Parses the modifications present in the sequence content and generates a list of
    modifications specific to the provided sequence type. Depending on whether the
    sequence type is protein or rna/dna, it creates either a `ResidueModification`
    or `NucleotideModification` object and collects these into a list.

    Parameters
    ----------
    seq_type : SequenceType
        The type of the sequence, which indicates whether the sequence contains
        protein or nucleotide-related data.
    seq_content : dict
        Dictionary containing sequence details. The key 'modifications' within this
        dictionary holds a list of modification records, each specifying the type of
        modification and its associated position.

    Returns
    -------
    list[Modification]
        A list of modifications generated based on the sequence type and content.
    """
    modifications = []
    for modification in seq_content["modifications"]:
        if seq_type == SequenceType.PROTEIN:
            tmp_mod = ResidueModification(
                mod_str=modification["ptmType"],
                mod_pos=modification["position"]
            )
        else:
            tmp_mod = NucleotideModification(
                mod_str=modification["modificationType"],
                mod_pos=modification["basePosition"]
            )
        modifications.append(tmp_mod)
    return modifications


def _parse_templates(
        seq_type: SequenceType,
        seq_content: dict
) -> list[Template]:
    """
    Parses sequence templates based on the provided sequence type and content.

    Parameters
    ----------
    seq_type : SequenceType
        The type of the sequence being processed.
    seq_content : dict
        A dictionary containing sequence template data. Includes details
        such as template definitions (strings or files) and indices
        for query and target alignment.

    Returns
    -------
    list[Template]
        A list of parsed `Template` objects constructed from the given
        sequence content.

    Raises
    ------
    AFTemplateError
        If templates are not supported for the given sequence type.
    ValueError
        If the type key of a template is invalid.
    """
    templates = []
    if seq_type != SequenceType.PROTEIN:
        raise AFTemplateError(
            f"Templates are not supported for {seq_type.value}."
        )
    for template in seq_content["templates"]:
        if TemplateType.STRING.value in template.keys():
            template_type = TemplateType.STRING
        elif TemplateType.FILE.value in template.keys():
            template_type = TemplateType.FILE
        else:
            raise ValueError("Invalid template structure.")
        tmp_template = Template(
            template_type=template_type,
            mmcif=template[template_type.value],
            qidx=template["queryIndices"],
            tidx=template["templateIndices"]
        )
        templates.append(tmp_template)
    return templates


def _parse_msa(
        seq_type: SequenceType,
        seq_content: dict
) -> MSA | None:
    """
    Parses multiple sequence alignments (MSA) based on the input sequence
    type and sequence content.

    This static method processes sequence content to construct a multiple
    sequence alignment (MSA) object if applicable. The method considers both
    paired and unpaired alignments provided either directly or via file paths.

    Parameters
    ----------
    seq_type : SequenceType
        The type of the sequence being processed.
    seq_content : dict
        A dictionary containing information about the sequence alignments.
        Possible keys include "pairedMsa", "pairedMsaPath", "unpairedMsa",
        and "unpairedMsaPath".

    Returns
    -------
    Optional[MSA]
        An MSA object if the input sequence content supports it, otherwise `None`.
    """
    msa = None

    paired = seq_content.get("pairedMsa") or seq_content.get("pairedMsaPath")
    unpaired = seq_content.get("unpairedMsa") or seq_content.get("unpairedMsaPath")

    paired_path = "pairedMsaPath" in seq_content.keys()
    unpaired_path = "unpairedMsaPath" in seq_content.keys()

    if seq_type == SequenceType.DNA and paired is not None:
        raise AFMSAError("Paired NSA is not supported for RNA sequences.")

    if paired is not None or unpaired is not None:
        msa = MSA(paired, unpaired, paired_path, unpaired_path)
    return msa


def _parse_sequence(seq_type: str, seq_content: dict) -> Sequence:
    """
    Parses a sequence from the given sequence type and content.

    This method processes a sequence from AlphaFold3 input data, including its type,
    modifications, templates, multiple sequence alignment (MSA).

    Parameters
    ----------
    seq_type : str
        The type of the sequence, which should correspond to a valid sequence type.
    seq_content : dict
        Dictionary containing the sequence entry data.

    Returns
    -------
    Sequence
        A `Sequence` object containing the parsed sequence data, along with its
        type, ID, modifications, templates, and MSA.
    """
    seq_id = _get_id(seq_content)
    seq_type = SequenceType(seq_type)

    modifications = None
    if "modifications" in seq_content.keys():
        modifications = _parse_modifcations(seq_type, seq_content)

    templates = None
    if "templates" in seq_content.keys():
        templates = _parse_templates(seq_type, seq_content)

    msa = None
    if seq_type != SequenceType.DNA:
        msa = _parse_msa(seq_type, seq_content)

    description = seq_content.get("description")

    return Sequence(
        seq_type=seq_type,
        seq_str=seq_content["sequence"],
        seq_id=seq_id,
        description=description,
        modifications=modifications,
        templates=templates,
        msa=msa,
    )


def read_json(filename: str, check: bool = True) -> InputFile:
    """
    Reads a JSON file and constructs an `InputFile` object by parsing the data.

    This method reads a given JSON file to extract information about the input
    system components, such as sequences, ligands, and bonds, and creates an
    `InputFile` object using the parsed data. It also performs some data consistency
    checks if specified.

    Parameters
    ----------
    filename : str
        Path to the input JSON file containing the system description.
    check : bool, optional
        Whether to perform a data consistency check after reading the input data.
        Default is True.

    Returns
    -------
    InputFile
        A fully constructed InputFile object based on the input JSON data.

    Raises
    ------
    AFMissingFieldError
        If an invalid bonded atom pair is detected in the JSON data, specifically
        when a pair does not contain exactly two elements or is malformed.
    """
    data = _read(filename)
    if check:
        _check_data(data)

    builder = InputBuilder()
    builder.set_name(data.get("name", ""))\
        .set_version(data.get("version", 1))\
        .set_dialect(data.get("dialect", "alphafold3"))\
        .set_seeds(data.get("modelSeeds", [1]))\

    for sequence in data["sequences"]:
        seq_type, seq_content = next(iter(sequence.items()))

        if seq_type == "ligand":
            ligand = _parse_ligand(seq_content)
            builder.add_ligand(ligand)
            continue

        next_sequence = _parse_sequence(seq_type, seq_content)
        builder.add_sequence(next_sequence)

    if "bondedAtomPairs" in data and data["bondedAtomPairs"] is not None:
        for pair in data["bondedAtomPairs"]:
            if len(pair) != 2:
                raise AFMissingFieldError(
                    f"Invalid bonded atom pair: {pair}"
                )
            atom1 = Atom(*pair[0])
            atom2 = Atom(*pair[1])
            bond = Bond(atom1, atom2)
            builder.add_bonded_atom_pair(bond)

    if "userCCD" in data:
        builder.set_user_ccd(data["userCCD"])

    return builder.build()
