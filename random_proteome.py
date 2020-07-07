"""Create random proteome from template"""

import argparse
import os
import random
import sys

from Bio import SeqIO


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
        List of proteins name in order of appearance in the fasta file.
    dict
        Dictionnary with protein names as keys and protein sequences as keys.
    """
    proteins = []
    proteome = {}
    with open(fasta_filename, "r") as proteome_file:
        for record in SeqIO.parse(proteome_file, "fasta"):
            name = record.id
            proteins.append(name)
            proteome[name] = str(record.seq)
    return proteins, proteome


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


def calculate_amino_acid_proportion(amino_acid_dct):
    """
    Calculate proportion of amino acids.

    Parameters
    ----------
    amino_acid_dct : dict
        Dictionnary with amino acid as key and count as value.

    Returns
    -------
    dict
        Dictionnary with amino acid as key and proportion as value.
    """
    amino_acid_prop = {}
    total = sum(amino_acid_dct.values())
    for amino_acid in amino_acid_dct:
        amino_acid_prop[amino_acid] = amino_acid_dct[amino_acid] / total
    return amino_acid_prop
    


def create_random_protein(length, amino_acid_lst, amino_acid_dct,
                          start_with_methionine=True,  
                          amino_acid_bias=None):
    """

    https://docs.python.org/3/library/random.html#random.choices
    """
    sequence = ""
    while len(sequence) != length:
        if (start_with_methionine == True 
            and sequence == "" 
            and amino_acid_dct["M"] > 0):
            amino_acid = "M"
        else:
            amino_acid = random.choices(amino_acid_lst, 
                                        weights=amino_acid_bias)[0]
        if amino_acid_dct[amino_acid] > 0:
            sequence += amino_acid
            amino_acid_dct[amino_acid] -= 1
    return sequence, amino_acid_dct
    

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
