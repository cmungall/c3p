"""
Classifies: CHEBI:75769 B vitamin
"""
from rdkit import Chem
from rdkit.Chem import AllChem


def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is likely a B vitamin based on its SMILES string.
    Checks for key structural features of B1, B2, B3, B5, B6, B7 and B9 vitamins.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a B vitamin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # B1 (Thiamine) - Thiazole and pyrimidine rings
    thiamine_pattern = Chem.MolFromSmarts('c1[n+]cscc1-CC-c1cncc(N)c1')
    if mol.HasSubstructMatch(thiamine_pattern):
        return True, "Contains thiamine (B1) substructure."

    # B2 (Riboflavin) - Isoalloxazine ring system
    riboflavin_pattern = Chem.MolFromSmarts('C1=CC=C2N(C[C@H](O)[C@H](O)[C@H](O)CO)C3=NC(=O)NC(=O)C3=NC2=C1')
    if mol.HasSubstructMatch(riboflavin_pattern):
          return True, "Contains riboflavin (B2) substructure."


    # B3 (Niacin/Nicotinamide) - Pyridine with carboxyl/amide
    nicotinamide_pattern = Chem.MolFromSmarts('c1ccncc1C(=O)[N]')
    nicotinic_acid_pattern = Chem.MolFromSmarts('c1ccncc1C(=O)O')
    if mol.HasSubstructMatch(nicotinamide_pattern) or mol.HasSubstructMatch(nicotinic_acid_pattern):
          return True, "Contains niacin/nicotinamide (B3) substructure."

    # B5 (Pantothenic acid) - Pantoic acid and beta-alanine
    pantothenic_pattern = Chem.MolFromSmarts('CC(C)(CO)[C@H](O)C(=O)NCCC(O)=O')
    if mol.HasSubstructMatch(pantothenic_pattern):
        return True, "Contains pantothenic acid (B5) substructure."

    # B6 (Pyridoxine, Pyridoxal, Pyridoxamine) - Pyridine with OH, CH2OH, CHO, and/or CH2NH2 groups
    pyridoxine_pattern = Chem.MolFromSmarts('C1=C(O)C(CO)=C(CO)C=N1')
    pyridoxal_pattern = Chem.MolFromSmarts('C1=C(O)C(CO)=C(C=O)C=N1')
    pyridoxamine_pattern = Chem.MolFromSmarts('C1=C(O)C(CO)=C(CN)C=N1')

    if mol.HasSubstructMatch(pyridoxine_pattern) or mol.HasSubstructMatch(pyridoxal_pattern) or mol.HasSubstructMatch(pyridoxamine_pattern):
        return True, "Contains pyridoxine/pyridoxal/pyridoxamine (B6) substructure."

    # B7 (Biotin) - Tetrahydrothiophene and imidazole
    biotin_pattern = Chem.MolFromSmarts('C1S[CH]2[CH]([CH2][C](N2)=O)NC(=O)N[C@H]3CCCC[C@H]3')
    if mol.HasSubstructMatch(biotin_pattern):
        return True, "Contains biotin (B7) substructure."

    # B9 (Folic acid/folate) - Pteridine with PABA and glutamic acid
    folic_acid_pattern = Chem.MolFromSmarts('c1nc2nc(N)[nH]c(=O)c2n1-CC-c1ccc(N)cc1-C(=O)NC(CCC(=O)O)C(=O)O')
    if mol.HasSubstructMatch(folic_acid_pattern):
        return True, "Contains folic acid (B9) substructure."



    return False, "Not recognized as a B vitamin based on substructure search."