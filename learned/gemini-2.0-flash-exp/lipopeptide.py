"""
Classifies: CHEBI:46895 lipopeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide contains a peptide (amino acid chain) and a lipid (long hydrocarbon chain).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Peptide backbone: Look for multiple amide bonds
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])N") # C=O-N
    amide_matches = mol.GetSubstructMatches(amide_pattern)

    if len(amide_matches) < 2: # Lipopeptides contain a chain at least 2 amino acids.
        return False, "Too few amide bonds, not a peptide"

    # 2. Lipid component: Look for long aliphatic chain 
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    
    if len(chain_matches) == 0:
        return False, "No long aliphatic chain found"

    # 3. Check for a linkage between lipid and peptide
    #  -Check for connectivity between chain and amide
    linked_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX3](=[OX1])N") #  C-C(=O)N 
    linked_chain_matches = mol.GetSubstructMatches(linked_chain_pattern)

    linked_O_pattern = Chem.MolFromSmarts("[CX4,CX3]~[OX2][CX3](=[OX1])") # C-O-C=O
    linked_O_matches = mol.GetSubstructMatches(linked_O_pattern)
    if not linked_chain_matches and not linked_O_matches:
        return False, "No connectivity between peptide and lipid"

    # 4. Molecular weight check. Lipopeptides are larger molecules, typically.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400: # Lipopeptides tend to be relatively large.
        return False, "Molecular weight too low for a lipopeptide"

    # Additional check - number of carbons and oxygens
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if carbon_count < 15:
         return False, "Too few carbon atoms for a lipopeptide"

    if oxygen_count < 3:
        return False, "Too few oxygen atoms for a lipopeptide"
    
    return True, "Contains a peptide chain linked to a lipid chain"