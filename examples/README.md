# af3cli Examples

This directory contains example files demonstrating the usage of `af3cli` through both CLI inputs and as a Python library. 

> [!NOTE]
> The CLI examples use associative arrays to improve clarity, therefore at least Bash version 4 is required.

## Example Directories and Purpose
The examples are organized into these main directories:
1. `cli/` - Demonstrates how to use `af3cli` via the command-line interface.
2. `python/` - Covers examples for using `af3cli` as a Python library in various workflows.
3. `jupyter/` - Provides interactive examples using Jupyter notebooks.
4. `data/` - Contains any additional files for sequences and ligands needed to run the examples

These directories contain examples illustrating the same core scenarios for consistency.

## Examples Overview
1. **Bonds**
   - Assigning a covalent bonds to a residue in a protein sequence utilizing a Ligand molecule in CCD format.
2. **DNA**
   - A simple example with double-stranded DNA (dsDNA).
3. **Modifications and PTMs**
   - Handling modifications of DNA/RNA and protein post-translational modifications (PTMs).
4. **Protein Systems with Ligands**
   - Setting up a protein system with a Ligand (no pre-defined interactions):
     - Using Ligands in **SMILES** format.
     - Using Ligands in **CCD** format.
5. **FASTA Reader**
   - Parsing and reading FASTA files using the built-in reader functions.

