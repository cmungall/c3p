"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: CHEBI:35368 fatty acid methyl ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is a fatty acid ester that is the carboxylic ester 
    obtained by the formal condensation of a fatty acid with methanol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for methyl ester group (-C(=O)O-C)
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, "Must have exactly one methyl ester group"

    # Look for long carbon chain (fatty acid part)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, "Missing fatty acid chain"

    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 7:
        return False, "Chain too short to be a fatty acid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 8:
        return False, "Too few carbons for fatty acid methyl ester"
    if o_count != 2:
        return False, "Must have exactly 2 oxygens (methyl ester)"

    # Check if the carbon chain is attached to the ester group
    ester_atom_idx = ester_matches[0][0]
    ester_atom = mol.GetAtomWithIdx(ester_atom_idx)
    for neighbor_atom in ester_atom.GetNeighbors():
        if neighbor_atom.GetAtomicNum() == 6:
            break
    else:
        return False, "Fatty acid chain not attached to ester group"

    return True, "Contains a methyl ester group attached to a fatty acid chain"