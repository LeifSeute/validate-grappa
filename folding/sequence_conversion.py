def three_to_one_letter(sequence):
    # Amino acid conversion dictionary
    aa_dict = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    }

    # Remove dashes and split the sequence into three-letter codes
    parts = sequence.split('-')
    
    # Convert to single-letter codes
    single_letter_sequence = ''.join(aa_dict[part] for part in parts)

    return single_letter_sequence

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='Convert a sequence from three-letter to one-letter code')
    parser.add_argument('--sequence', help='The sequence to convert', type=str)
    args = parser.parse_args()
    sequence = args.sequence

    print(f"Converting sequence {sequence} to one-letter code")
    print(three_to_one_letter(sequence))


# paper trp cage: ASP-ALA-TYR-ALA-GLN-TRP-LEU-ALA-ASP-GLY-GLY-PRO-SER-SER-GLY-ARG-PRO-PRO-PRO-SER
# paper villin: LEU-SER-ASP-GLU-ASP-PHE-LYS-ALA-VAL-PHE-GLY-MET-THR-ARG-SER-ALA-PHE-ALA-ASN-LEU-PRO-LEU-TRP--NLE-GLN-GLN-HIS-LEU-NLE-LYS-GLU-LYS-GLY-LEU-PHE
# paper a3d: MET-GLY-SER-TRP-ALA-GLU-PHE-LYS-GLN-ARG-LEU-ALA-ALA-ILE-LYS-THR-ARG-LEU-GLN-ALA-LEU-GLY-GLY-SER-GLU-ALA-GLU-LEU-ALA-ALA-PHE-GLU-LYS-GLU-ILE-ALA-ALA-PHE-GLU-SER-GLU-LEU-GLN-ALA-TYR-LYS-GLY-LYS-GLY-ASN-PRO-GLU-VAL-GLU-ALA-LEU-ARG-LYS-GLU-ALA-ALA-ALA-ILE-ARG-ASP-GLU-LEU-GLN-ALA-TYR-ARG-HIS-ASN
    
# paper protein g: ASP-THR-TYR-LYS-LEU-VAL-ILE-VAL-LEU-ASN-GLY-THR-THR-PHE-THR-TYR-THR-THR-GLU-ALA-VAL-ASP-ALA-ALA-THR-ALA-GLU-LYS-VAL-PHE-LYS-GLN-TYR-ALA-ASN-ASP-ALA-GLY-VAL-ASP-GLY-GLU-TRP-THR-TYR-ASP-ALA-ALA-THR-LYS-THR-PHE-THR-VAL-THR-GLU
    
# ww domain: GLY-SER-LYS-LEU-PRO-PRO-GLY-TRP-GLU-LYS-ARG-MET-SER-ARG-ASP-GLY-ARG-VAL-TYR-TYR-PHE-ASN-HIS-ILE-THR-GLY-THR-THR-GLN-PHE-GLU-ARG-PRO-SER-GLY

# ntl9: MET-LYS-VAL-ILE-PHE-LEU-LYS-ASP-VAL-LYS-GLY-MET-GLY-LYS-LYS-GLY-GLU-ILE-LYS-ASN-VAL-ALA-ASP-GLY-TYR-ALA-ASN-ASN-PHE-LEU-PHE-LYS-GLN-GLY-LEU-ALA-ILE-GLU-ALA (CONH2)