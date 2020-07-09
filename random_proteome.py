"""Create random proteome from template"""

import argparse
import os
import random
import sys

from Bio import SeqIO
import pandas as pd

random.seed(42)
AMINO_ACID_LST = list("ACDEFGHIKLMNPQRSTVWY")



def cli():
    """Get options from commande line inteface.

    Parameters
    ----------

    Returns
    -------

    """
    parser = argparse.ArgumentParser(
        description="Create random proteome from template."
    )
    # Define arguments.
    parser.add_argument("-t", "--template", action="store", dest="template",
                        required=True,
                        help="Name of template proteome file.")
    # Parse arguments.
    inputs = parser.parse_args()
    # Verify file exits.
    if not os.path.isfile(inputs.template):
        sys.exit(f"{inputs.template} is not a valid file. Bye.")
    return inputs


def read_template_proteome(fasta_filename):
    """Read template proteome from file.

    Parameters
    ----------
    fasta_file : str
        Name of FASTA file containing template proteome.

    Returns
    -------
    list
        List of protein names in order of appearance in the fasta file.
    list
        List of protein sequences in order of appearance in the fasta file.
    """
    protein_names = []
    protein_sequences = []
    with open(fasta_filename, "r") as proteome_file:
        for record in SeqIO.parse(proteome_file, "fasta"):
            protein_names.append(record.id)
            protein_sequences.append(str(record.seq))
    return protein_names, protein_sequences


def count_amino_acids(sequence, amino_acid_lst):
    """Count amino acids.
    """
    count_dct = {}
    for amino_acid in amino_acid_lst:
        count_dct[amino_acid] = sequence.count(amino_acid)
    return count_dct


def add_dictionnaries(dict_1, dict_2):
    """Add to dictionnaries.

    Assumes both dictionnaries have the same keys.

    Parameters
    ----------
    dict_1 : dict
        First dictionnary.
    dict_2 : dict
        Second dictionnary.
    
    Returns
    -------
    dict
        Sum of both dictionnaries.
    """
    assert sorted(sorted(dict_1.keys())) == sorted(sorted(dict_2.keys()))
    dict_sum = {}
    for key in dict_1:
        dict_sum[key] = dict_1[key] + dict_2[key]
    return dict_sum


def get_amino_acid_proportion(sequence, amino_acid_lst):
    """
    Calculate proportion of amino acids.

    Parameters
    ----------
    amino_acid_lst : list
        List of amino acids to consider.

    Returns
    -------
    dict
        Dictionnary with amino acid as key and proportion as value.
    """
    amino_acid_prop = {}
    for amino_acid in amino_acid_lst:
        amino_acid_prop[amino_acid] = sequence.count(amino_acid) / len(sequence)
    return amino_acid_prop


def shuffle_sequence(sequence):
    """Shuffle sequence.
    
    Parameters
    ----------
    sequence : str
        Sequence to shuffle.
    
    Returns
    -------
    str
        Shuffled sequence.
    """
    shuffled_sequence = list(sequence)
    random.shuffle(shuffled_sequence)
    return "".join(shuffled_sequence)


def create_random_protein_from_proteome(proteome_sequence, length):
    """Create random protein from a random proteome.
    
    Parameters
    ----------
    proteome_sequence : str
        Complete sequence of the entire proteome.
    length : int
        Lenght of the target protein.
    
    Returns
    -------
    str
        Produced protein sequence.
    str
        Remaining proteome sequence.
    """
    protein_sequence = proteome_sequence[0:length]
    proteome_sequence = proteome_sequence[length:]
    return protein_sequence, proteome_sequence


def create_random_protein_from_distribution(length, amino_acid_lst, amino_acid_distribution=None):
    """Create random protein from distribution.

    https://docs.python.org/3/library/random.html#random.choices
    
    Parameters
    ----------
    lenght : int
        Lenght of the target protein.
    amino_acid_lst : list
        List of amino acids to consider.
    amino_acid_distribution : list
        List of amino acids distribution.
        
    Returns
    -------
    str
        Produced protein sequence.
    """
    sequence = random.choices(amino_acid_lst, 
                              weights=amino_acid_distribution,
                              k=length)
    return "".join(sequence)


def write_fasta(protein_lst, fasta_filename, width=60):
    """Write sequence as fasta file.

    Parameters
    ----------
    protein_lst : list
        List of proteins.
    fasta_filename : str
        Name of output fasta file.
    width : int, optional
        Width of sequence line in fasta format.
    """
    with open(fasta_filename, "w") as fasta_file:
        for prot_idx, prot_seq in enumerate(protein_lst):
            size = len(prot_seq)
            # Protein index starts at 1 in the description line.
            fasta_file.write(f">protein: {prot_idx+1} | size: {size}\n")
            for i in range(0, len(prot_seq), width):
                fasta_file.write(f"{prot_seq[i: i+width]}\n")


def write_distribution(protein_lst, amino_acid_lst, filename, ref_distribution=None):
    """Write amino acid distribution of proteins.

    Parameters
    ----------
    protein_lst : list
        List of proteins.
    amino_acid_lst : list
        List of reference amino acids.
    filename : str
        Name of output file in .tsv format.
    ref_distribution : dict, optional
        Dictionnary with amino acid ditribution of the reference proteome.
    """
    df  = pd.DataFrame(columns=amino_acid_lst)
    if ref_distribution:
        df = df.append(pd.Series(ref_distribution, name=0))
    for prot_index, prot_sequence in enumerate(protein_lst):
        proportion = get_amino_acid_proportion(prot_sequence, amino_acid_lst)
        df = df.append(pd.Series(proportion, name=prot_index+1))
    df.to_csv(filename, sep="\t", index_label="prot_id")

if "__name__" == "__main__":
    inputs = cli()
    print(f"Template proteome is:{inputs.template}")
    PROTEIN_NAME_LST, PROTEOME_DCT = read_template_proteome(inputs.template)
    
    sys.exit(1)
    amino_acid_proteome_count = dict(zip(list(AMINO_ACID_STR), [0]*20))
    # Count amino acids in template proteome.
    for protein_name in PROTEIN_NAME_LST:
        amino_acid_protein_count = count_amino_acids(
                                        PROTEOME_DCT[protein_name],
                                        AMINO_ACID_LST)
        amino_acid_proteome_count = add_dictionnaries(
                                        amino_acid_proteome_count, 
                                        amino_acid_protein_count)

    # Build random proteome.
