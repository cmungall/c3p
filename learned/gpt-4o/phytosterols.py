"""
Classifies: CHEBI:26125 phytosterols
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols similar to cholesterol with variations mainly in carbon side chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the sterol core pattern (four-ring structure typical of sterols)
    sterol_pattern = Chem.MolFromSmarts('C1CCC2C3CCC4CCCC(C4)C3C2C1')
    if not mol.HasSubstructMatch(sterol_pattern):
        return False, "No sterol core found"

    # Identify variations in the side chain (e.g., hydrocarbon chains with or without double bonds)
    # This is more flexible compared to a precise match for each phytosterol.
    side_chain_pattern = Chem.MolFromSmarts('CCCC=C|CC=C|CCC')
    side_chain_matches = mol.GetSubstructMatches(side_chain_pattern)
    if len(side_chain_matches) < 1:
        return False, "Missing characteristic side chain extensions"

    # Check for the number of carbons and presence of alcohol group
    has_alcohol = any(atom.GetAtomicNum() == 8 and atom.GetTotalDegree() == 1 for atom in mol.GetAtoms())
    if not has_alcohol:
        return False, "No hydroxyl group attached to the sterol core"

    # Check for a double bond in the side chain for variation
    has_double_bond = any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE for bond in mol.GetBonds())
    if has_double_bond:
        return True, "Contains sterol core with side chain variations and double bonds"

    return True, "Contains sterol core with characteristic side chains"

# Example usage:
# smiles = "SMILES_OF_A_CHEMICAL"
# print(is_phytosterols(smiles))