"""
Example:
tRNA-mimicking RNA from TYM Virus
    - Protein Data Bank (PDB) Entry: 4P5J
    - Reference: https://www.rcsb.org/structure/4P5J

This script provides a Python implementation as an alternative to the af3cli input bash script.
"""

import pprint
from pathlib import Path
from af3cli import InputBuilder, RNASequence
from af3cli.sequence import read_fasta

# Define constants
FILENAME = "example_FASTA_reader_python.json"
JOB_NAME = "example_FASTA_reader_py_job"

INPUT_FASTA_FILEPATH = str(
    Path(__file__).resolve().parent / "../data/rcsb_pdb_4P5J.fasta"
)

# Build input configuration for the job
input_builder = InputBuilder()
input_builder.set_name(JOB_NAME)

# Create RNA sequence object
for name, fasta_sequence_string in read_fasta(INPUT_FASTA_FILEPATH):
    sequence = RNASequence(
        seq_str=fasta_sequence_string,
        description=name,
    )
    input_builder.add_sequence(sequence)

internal_input = input_builder.build()

# Uncomment following line to generate output as JSON file
#internal_input.write(FILENAME)

print_json_via_debug = pprint.PrettyPrinter(indent=4)
print_json_via_debug.pprint(internal_input.to_dict())
