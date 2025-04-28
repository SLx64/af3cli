# af3cli
A command-line interface and Python library for generating [AlphaFold3](https://github.com/google-deepmind/alphafold3) input files.

## Installation

We recommend using [uv](https://github.com/astral-sh/uv) to manage your installation.

```shell
git clone https://github.com/SLx64/af3cli.git
cd af3cli
uv sync --locked
```

This automatically creates a virtual environment `.venv` in the project folder and installs all dependencies. If you do not need the optional dependencies for reading SDF ([RDKit](https://github.com/rdkit/rdkit)) or FASTA files ([Biopython](https://github.com/biopython/biopython)), the installation can be prevented with `--no-group features`.

The latest release is also available via [PyPI](https://pypi.org/project/af3cli/) and can be installed with `pip`. The optional dependencies must be specified here explicitly.

```shell
pip install af3cli [biopython rdkit]
```

## Citation

Please cite the following article when using af3cli in your published work:

> Döpner, P.; Kemnitz, S.; Doerr, M.; Schulig, L.; af3cli: Streamlining AlphaFold3 Input Preparation. Journal of Chemical Information and Modeling **2025**, *65* (8), 3886-3891. DOI: [10.1021/acs.jcim.5c00276](https://doi.org/10.1021/acs.jcim.5c00276) 

## Basic Usage

The generation of AlphaFold3 input files can be done either with the standalone CLI tool or for more advanced tasks by using the library in Python scripts. [Python Fire](https://github.com/google/python-fire) is used to implement the CLI application.

For a detailed overview of all available JSON fields check the [AlphaFold3 input documentation](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md).

> [!WARNING]
> In most cases, checks are only carried out to ensure a correct structure of the input file, but not whether the inputs themselves are valid.

You can display the help and overview of the available CLI commands with the following statement.

```shell
af3cli -- --help
```

All commands and sub-commands are separated with `-` to enable chaining, allowing e.g. several sequences, ligands, bonds etc. to be added.

```shell
af3cli toplevel sub [...] - sub [...] \
    - toplevel sub [...]
```

You can use the `debug` command to display the final file without writing it.

```shell
af3cli debug --show - [...]
```

### Config Parameters

The `config` command is used to manage basic settings, such as the file name of the JSON file to be written, the name of the job or the respective version.

```shell
af3cli config -f "filename.json" -j "jobname" -v 2
```

The library provides an `InputBuilder` that allows new jobs to be created very comfortably step by step.

```python
from af3cli import InputBuilder

builder = InputBuilder()
builder.set_name("jobname")
builder.set_version(2)
# builder.set_dialect("alphafold3") # default

input_file = builder.build()
input_file.write("filename.json")
```

You can also initialize the `InputBuilder` with an existing `InputFile` object in order to add further sequences or ligands or to change settings.

### Random Seeds

It is required that at least one random seed is specified. The default value is therefore 1. Otherwise, you can either specify a number of values to generate a list of random seeds or pass a list of integers yourself.

```shell
af3cli seeds -n 10 - ...
# generates 10 random numbers

af3cli seeds -v "1,2,3"
# "(1,2,3,...)" or "[1,2,3,...]" are also valid
```

Python:
```python
builder.set_seeds([1, 2, 3])
```

### Sequences

Adding sequences works basically the same for all three available types, but not all JSON fields are supported for each type. The corresponding subcommands therefore differ in some cases.

```shell
af3cli [...] \
    - protein add "MVKLAGST" \ # positional argument
    - protein add --sequence "AAQAA" \
    - dna add --sequence "AATTTTCC" \
    - rna add --sequence "UUUGGCCGG"
```

A check is performed to ensure that the sequence characters match the respective type in the CLI application or in the Python library itself, when the `Sequence` object is converted into a dictionary. When using the sequence base class, the corresponding `SequenceType` must be specified. We therefore provide derived classes for the respective sequence types to simplify use.

```python
from af3cli import ProteinSequence # DNASequence, RNASequence

# DNASequence / RNASequence
protein_seq = ProteinSequence(
    "MVKLAGST...",
    #[...]
)

builder.add_sequence(protein_seq)
```

For DNA sequences, it is possible to generate the reverse complementary strand. The associated data, such as manually specified IDs or modifications, are not included. The CLI tool will generate appropriate warnings.

```shell
af3cli [...] \
    - dna add --sequence "AATTTTCC" --complement
```

Python:

```python
from af3cli import DNASequence


dna_seq = DNASequence(
    "AATTTTCC...",
    #[...]
)
rc_dna_seq = dna_seq.reverse_complement()

builder.add_sequence(dna_seq)
builder.add_sequence(rc_dna_seq)
```

If modifications or manually defined IDs are required, the complementary sequence must be created separately.

#### FASTA Files

As it is often not very practical to add many or particularly long sequences via the CLI, it is possible to read the respective sequence from a FASTA file. To use this feature, [Biopython](https://github.com/biopython/biopython) must be installed as an optional dependency.

```shell
af3cli protein add --sequence <filename> --fasta
```

Each sequence command expects exactly one single sequence. Otherwise it is not possible to add additional fields, such as modifications or templates. However, it is still possible to read several sequences from a FASTA file if the additional features are not required. For even more advanced tasks, the Python API must otherwise be used.

```shell
af3cli [...] - fasta [--filename] <filename>
```

The respective sequence type is automatically detected, which is not possible in rare cases. If this is the case, the sequence is ignored and a warning is issued. It is, therefore, advisable to add all sequences whose type cannot be clearly identified separately via the sequence commands.

There are also two ways of doing this when using Python. The `fasta2seq` function can be used to obtain a generator that automatically creates `Sequence` objects and the `read_fasta` function is used to create a generator that returns the plain FASTA IDs and sequences from the FASTA file as a string.

```python
from af3cli.sequence import fasta2seq, read_fasta

for seq in fasta2seq(filename):
    ...
    # do something with the Sequence object

for fasta_id, seq_str in read_fasta(filename):
    ...
    # create your own Sequence objects
```

#### Modifications

By applying the `modification` subcommand, any number of modifications can be added to the sequences with the respective CCD identifier and position. The different fields in the JSON file are automatically inserted correctly based on the sequence type.

```shell
# as positional arguments
af3cli [...] protein [...] - modification "SEP" 5
# or with explicit argument names
af3cli [...] dna [...] - modification --mod "6OG" --pos 1
```

When using the Python API, you have to explicitly define what kind of modifications you would like to add, since the resulting JSON fields are different for protein (`ResidueModification`) or nucleotide sequences (`NucleotideModification`). Please note that checks are performed to verify the modification types when the `Sequence` object is converted to a dictionary.

```python
from af3cli import ProteinSequence, ResidueModification
                   # NucleotideModification

rmod = ResidueModification("SEP", 5)

protein_seq = ProteinSequence(
    "<SEQUENCE>",
    # [...]
    modifications=[rmod]
)

# it is possible to add more modifications later
protein_seq.modifications.append(rmod)
```

#### Structural Templates

For protein sequences, it is possible to specify multiple structural templates in mmCIF format as a string or path. Since it is completely impractical to use strings via the CLI tool, the file must be submitted as plain text and is then read in its entirety as a string. 

```shell
# read the file as string with the '--read' flag
af3cli [...] protein [...] - template [--mmcif] <filename> --read
# keep relative/absolute path
af3cli [...] protein [...] - template [--mmcif] <filename>

# specify query and template indices as list of integers
# "1,2,3,..." | "(1,2,3,...)" | "[1,2,3,...]" are valid
af3cli [...] protein [...] \
    - template [--mmcif] <filename> -q "..." -t "..."
```

As it makes no difference to Python whether the string contains a path to a file or the file content, all you need to do is specify the template type. The file must then be read manually beforehand if a string is desired in the JSON file.

```python
from af3cli import Template, TemplateType, ProteinSequence

# TemplateType.FILE for relative/absolute path
t = Template(
    TemplateType.STRING,
    "mmCIF content",
    qidx=[], tidx=[]
)

protein_seq = ProteinSequence(
    "<SEQUENCE>",
    # [...]
    templates=[t]
)

# it is possible to add more templates later
protein_seq.templates.append(t)
```

#### Multiple Sequence Alignment

Please refer to the [AlphaFold3 input documentation](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md#multiple-sequence-alignment) on how to specify the MSA section for protein and RNA sequences.

The A3M-formatted content can be specified either as a path or as a string (mutally exclusive).

```shell
af3cli [...] protein [...] msa --paired ... --unpaired ...
af3cli [...] protein [...] msa --pairedpath ... --unpairedpath ...
```

In the case of the Python API, you must specify whether the respective string is a path.

```python
from af3cli import MSA

msa = MSA(
    paired="...", unpaired="...",
    paired_is_path=True, unpaired_is_path=True,
)

protein_seq = ProteinSequence(
    "<SEQUENCE>",
    # [...]
    msa=msa
)

# alternative
protein_seq._msa = msa
```

### Ligands and Ions

The ligands are treated in a generally similar way to the sequences and can be defined either as SMILES or with a corresponding CCD identifier. SDF files can also be read and converted to SMILES via an optional [RDKit](https://github.com/rdkit/rdkit) dependency. If there are multiple entries in the SDF, they are added as individual ligands. Ions are simply treated as ligands in AlphaFold3.

```shell
af3cli [...] \
    - ligand add --smiles "CCC" \
    # providing a list of CCD codes is also supported
    - ligand add --ccd "MG" \
    - ligand add --sdf ligands.sdf
```

In Python, either the parent class `Ligand` together with the respective `LigandType` or alternatively the corresponding child classes `CCDLigand` or `SMILigand` can be used to add new ligands.

```python
from af3cli import Ligand, LigandType, SMILigand
from af3cli.ligand import sdf2smiles

ligand = Ligand(
    LigandType.SMILES,
    "CCC",
    #[...]
)

# using SMILigand
# ligand = SMILigand("CCC")

builder.add_ligand(ligand)

# ...
for smi in sdf2smiles("ligands.sdf"):
    builder.add_ligand(
        Ligand(LigandType.SMILES, smi)
    )
```

### Custom CCD

Please refer to the [AlphaFold3 input documentation](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md#user-provided-ccd-format) on how to generate valid CCD mmCIF files.

The entire file content is stored as a string in the JSON file and is only stored in a variable here. A plain text file must, therefore, simply be specified for the CLI.

```shell
af3cli [...] ccd [--filename] <filename>
```

In Python, you then have to read the file yourself.

```python
builder.set_user_ccd(filecontent)
```

### Bonds

The bonded atom pairs are defined in the JSON file as a list of lists, each of which contains the Entity ID, the Residue ID and the atom name. To make it as easy as possible to add new bonds, a string format is used, which is then translated into the correct format.

```shell
# E: Entity ID; R: Residue ID N: atom name
af3cli [...] bond [--add] "E:R:N-E:R:N"

# example
af3cli [...] bond [--add] "A:1:C-B:1:O"
```

Although the sequences should be numbered in the order in which they were added, it is advisable to manually assign a sequence ID to the respective entities for the bonds (see below).

Python:

```python
from af3cli import Bond

bond = Bond.from_string("A:1:C-B:1:O")

builder.add_bonded_atom_pair(bond)
```

You can also use the `Atom` class to initialize new atoms and create a `Bond` object from any two atoms.

```python
from af3cli import Bond, Atom

atom_1 = Atom("A", 1, "C")
atom_2 = Atom("B", 1, "O")
bond = Bond(atom_1,atom_2)

builder.add_bonded_atom_pair(bond)
```

### Sequence ID Handling

The IDs for sequences, ligands, and ions are normally assigned automatically and should only be specified manually if it is really necessary, as ID clashes may occur. An `IDRegister` object keeps track of the sequences used and, if necessary, skips IDs that have already been registered.

One case where it is necessary to specify the IDs manually is for bonds between different entries, as the chain ID must be specified for the bonded atom pairs (see above).

```shell
# "A,B,..." | "(A,B,...)" | "[A,B,...]" are valid
af3cli [...] protein add [...] -i "A,B"
```

If you only want to calculate homomultimers without specifying an explicit ID, you can also specify a number.

```shell
af3cli [...] protein add [...] -n 2

# works for all sequence types and ligands/ions
af3cli [...] \
    - protein add [...] -n 5 \
    - ligand add [...] -n 5
```

You can also specify IDs or a number in connection with an SDF file, whereby it should be noted that the number of manually specified IDs must correspond to the number of ligands in the SDF file. If a number is specified, all entries in the SDF are then multiplied by this number. Since all entries in the SDF file are converted to a SMILES string, manual ID assignment is not needed as bonded atoms are only available for CCD entries.

In Python, the number or explicit IDs can be specified when initializing `Ligand` or `Sequence` objects. If both parameters are specified, IDs are prioritized. If the number of IDs assigned is therefore greater than the number, the latter is overwritten. In the opposite case, missing IDs are populated automatically. This also facilitates subsequent changes to the number.

The registration or automatic assignment of IDs only takes place in connection with an `InputFile` object and is carried out when the file is converted into a dictionary (e.g. when the file is written). If the IDs of a sequence are changed after they have already been registered, the `IDRegister` must be reset.

```python
ligand = SMILigand(
    "CCC",
    seq_id=["A", "B"],
    num=2
)
```

### Merging Files

Occasionally, it can be helpful to create a base file of your system and prepare subsequent AlphaFold3 jobs by merging existing files with new entries. The `merge` command is chainable, allowing to combine several files. However, this should be done with caution if certain IDs, bonds, or seeds are important.

```shell
af3cli [...] merge [--filename] <filename>

# add new sequences
af3cli [...] merge [--filename] <filename> \
    - protein add "MVKLAGST..." \
    - ligand add --ccd "MG"

# keep IDs
af3cli [...] merge [--filename] <filename> --noreset

# override/merge special entries
af3cli [...] merge [--filename] <filename> \
    # override user-specified CCD data
    --userccd \
    # merge bonded atoms data
    --bonds \
    # merge seeds (removes duplicates)
    --seeds
```

Python:

```python
from af3cli import InputFile

input_file = InputFile()
other_input_file = InputFile.read("filename")

input_file.merge(other_input_file)

# with additional parameters
input_file.merge(
    other_input_file,
    reset=True,
    seeds=True,
    bonded_atoms=False,
    userccd=False
)
```
