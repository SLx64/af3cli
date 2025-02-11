# utility methods for convienence
import logging
from typing import Self

from af3cli.builder import InputBuilder
from af3cli.ligand import CCDLigand, SMILigand
from af3cli.sequence import SequenceType, Sequence, read_fasta, identify_sequence_type, ProteinSequence, DNASequence, RNASequence, sanitize_sequence_name
from af3cli.exception import AFSequenceTypeError, AFSequenceError

logger = logging.getLogger(__name__)

class SimpleBuilder(InputBuilder):
    """
    The SimpleBuilder class provides utility methods to 
    """

    def add_ligand_ccd(self, 
                   ccd_code: str = None,
                   seq_id: list[str] | None = None,
                   num: int = 1) -> Self:
        """
        Convenience method to add a single, CCD encoded ligand to the collection of 
        ligands in the `InputFile` instance.

        Parameters
        ----------
        ccd_code : str 
            The CCD (Chemical Component Dictionary) code of the ligand to be added. 
        seq_id : list[str], optional
            The sequence ID of the ligand to be added.
        num : int, optional
            The number of ligands to add. Defaults to 1.

        Returns
        -------
        Self
            Returns the current instance of the object to allow method chaining.
        """
        
        ligand = CCDLigand([ccd_code], seq_id=seq_id, num=num)
        self._afinput.ligands.append(ligand)
        return self
    
    def add_ligand_smiles(self, 
                   smiles_str: str | None = None,
                   seq_id: list[str] | None = None,
                   num: int = 1) -> Self:
        """
        Convenience method to add a single, SMILES encoded ligand to the collection of
        ligands in the `InputFile` instance.

        Parameters
        ----------
        smiles_str : str
            The SMILES string of the ligand to be added.
        seq_id : list[str], optional
            The sequence ID of the ligand to be added.
        num : int, optional
            The number of ligands to add. Defaults to 1.

        Returns
        -------
        Self
            Returns the current instance of the object to allow method chaining.
        """
       
        ligand = SMILigand([smiles_str], seq_id=seq_id, num=num)
        self._afinput.ligands.append(ligand)
        return self
    
    def add_sequence_fasta(self, fasta_filename: str = None,
                     num: int = 1,
                     add_num_to_all: bool = False,
                     sequence_type : SequenceType = None) -> Self:
        """
        Adds a sequence to the list of sequences in the `InputFile` instance.

        Parameters
        ----------
        sequence : Sequence
            The sequence to be added to the list of sequences.
        fasta_filename : str, optional
            The path to the FASTA file containing the sequence to be added.
        num : int, optional
            The number of sequences to add. Defaults to 1.
        add_num_to_all : bool, optional
            A flag to indicate whether the sequence number should be added to all sequences in the FASTA file.
            This is a safety feature to prevent the user from accidentally unintentionally adding the same number to all sequences.

        Returns
        -------
        Self
            Returns the current instance of the object to allow method chaining.
        """

        seq_tuples = [ seq_tuple for seq_tuple in read_fasta(fasta_filename)]

        if len(seq_tuples) > 1 and num > 1 and not add_num_to_all:
            logger.warning(f"Multiple sequences found in the FASTA file {fasta_filename}.")
            logger.warning(f"Please set add_num_to_all=True if you want to add the same number to all sequences.")
            raise AFSequenceError("Multiple sequences found in the FASTA file. Please set add_num_to_all=True if you want to add the same number to all sequences.")

        for seq_tuple in seq_tuples:
            if sequence_type is None:
                seq_type = identify_sequence_type(seq_tuple[1])
            if seq_type is None:
                raise AFSequenceTypeError("Could not automatically determine the sequence type., please specify it.")
            
            match seq_type:
                case SequenceType.PROTEIN:
                    sequence = ProteinSequence(seq_str=seq_tuple[1], seq_name=sanitize_sequence_name(seq_tuple[0]), num=num)
                case SequenceType.DNA:
                    sequence = DNASequence(seq_str=seq_tuple[1], seq_name=sanitize_sequence_name(seq_tuple[0]), num=num)
                case SequenceType.RNA:
                    sequence = RNASequence(seq_str=seq_tuple[1], seq_name=sanitize_sequence_name(seq_tuple[0]), num=num)
                
            self._afinput.sequences.append(sequence)
        return self

    
    def add_sequence_str(self, seq_str: str = None,
                         seq_type: SequenceType = None,
                     num: int = 1 ) -> Self:
        """
        Adds a sequence to the list of sequences in the `InputFile` instance.

        Parameters
        ----------
        seq_str : str, optional
            The string representation of the sequence to be added.
        num : int, optional
            The number of sequences to add. Defaults to 1.

        Returns
        -------
        Self
            Returns the current instance of the object to allow method chaining.
        """
        if seq_type is None:
            seq_type = identify_sequence_type(seq_str)
        if seq_type is None:
            raise AFSequenceTypeError("Could not automatically determine the sequence type., please specify it.")
        sequence = Sequence(seq_str=seq_str, seq_type=seq_type, num=num)
      
        self._afinput.sequences.append(sequence)
        return self

