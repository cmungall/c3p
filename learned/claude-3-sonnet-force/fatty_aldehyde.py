"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: CHEBI:35748 fatty aldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is an aldehyde formally arising from reduction of the carboxylic acid group
    of its corresponding fatty acid, having a carbonyl group at one end of the carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for aldehyde functional group (-C=O)
    aldehyde_pattern = Chem.MolFromSmarts("[CH3][CX3](=[OX1])")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde group found"
    
    # Look for long carbon chain
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if not chain_matches:
        return False, "No long carbon chain found"
    
    # Check if aldehyde is at the end of the chain
    aldehyde_atom = mol.GetSubstructMatch(aldehyde_pattern)[0]
    end_atom = None
    for match in chain_matches:
        if aldehyde_atom in match:
            end_atom = match[0] if match[0] != aldehyde_atom else match[-1]
            break
    if end_atom is None or mol.GetAtomWithIdx(end_atom).GetDegree() > 1:
        return False, "Aldehyde group not at the end of the chain"
    
    # Check for unsaturations in the carbon chain
    unsaturated = any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in mol.GetBonds())
    
    return True, "Contains aldehyde group at the end of a long carbon chain" + (", unsaturated" if unsaturated else "")