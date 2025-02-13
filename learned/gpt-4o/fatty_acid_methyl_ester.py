"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is characterized by having a methyl ester subgroup and generally
    a sufficient hydrophobic carbon chain.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")

    # Check for methyl ester group: C(=O)OC and ensure it's the terminating group
    methyl_ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    if not mol.HasSubstructMatch(methyl_ester_pattern):
        return False, "No methyl ester group found"
    
    # Look for long carbon chain characteristic of fatty acids
    # A more defined pattern for alkyl chains might be useful
    long_chain_pattern = Chem.MolFromSmarts("C(=O)O[C!H1]")  # CH3 at the end more restrictive pattern
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No sufficient hydrophobic tail detected"

    # Check the number of carbon atoms, typically FAMEs have more significant aliphatic carbon presence
    # Exact criteria might need fine-tuning based on molecular size/complexity judgments
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 8:
        return False, f"Insufficient carbon content for typical FAME; found {carbon_count} carbons"
    
    # Assess bond rotation to verify the chain structure
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, f"Insufficient rotatable bonds, indicating limited chain length or flexibility"

    return True, "Contains methyl ester group with an appropriate hydrophobic tail"