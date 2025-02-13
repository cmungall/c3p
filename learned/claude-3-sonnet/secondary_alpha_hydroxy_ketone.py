"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: CHEBI:33762 secondary alpha-hydroxy ketone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone has a carbonyl group and a hydroxy group linked
    by a carbon bearing one hydrogen and one organyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find carbonyl and hydroxy groups
    carbonyl_smarts = "[C=O]"
    hydroxy_smarts = "[OX2H]"
    carbonyl_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if mol.GetAtomWithIdx(atom.GetIdx()).HasQuery(Chem.MolFromSmarts(carbonyl_smarts))]
    hydroxy_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if mol.GetAtomWithIdx(atom.GetIdx()).HasQuery(Chem.MolFromSmarts(hydroxy_smarts))]
    
    # Check for at least one carbonyl and hydroxy group
    if not carbonyl_atoms or not hydroxy_atoms:
        return False, "Missing carbonyl or hydroxy group"
    
    # Find alpha carbon between carbonyl and hydroxy group
    alpha_carbons = []
    for c_idx in carbonyl_atoms:
        for h_idx in hydroxy_atoms:
            alpha_carbon = Chem.FindAtomEnvironmentOfRadiusN(mol, 2, c_idx, h_idx)
            if alpha_carbon:
                alpha_carbons.append(alpha_carbon[0].GetIdx())
    
    # Check if alpha carbon is secondary (has one H and one organyl group)
    for alpha_idx in set(alpha_carbons):
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        if alpha_atom.GetSymbol() == "C" and alpha_atom.GetTotalNumHs() == 1:
            if len(alpha_atom.GetNeighbors()) == 3:
                return True, "Contains a secondary alpha-hydroxy ketone moiety"
    
    return False, "No secondary alpha-hydroxy ketone moiety found"