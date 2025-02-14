"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: phospho sugar
"""
from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for phosphate group (considering various protonation states)
    phosphate_smarts = Chem.MolFromSmarts('OP(=O)(O)O')  # Phosphate group

    # Define SMARTS pattern for monosaccharide
    # This pattern looks for a sugar ring (5 or 6 membered ring with oxygens and carbons)
    sugar_smarts = Chem.MolFromSmarts('[#6,#8]1-[#6,#8]-[#6,#8]-[#6,#8]-[#6,#8]-1')

    # Find matches for phosphate group
    phosphate_matches = mol.GetSubstructMatches(phosphate_smarts)
    if not phosphate_matches:
        return False, "No phosphate group found"

    # Find matches for sugar moiety
    sugar_matches = mol.GetSubstructMatches(sugar_smarts)
    if not sugar_matches:
        return False, "No monosaccharide (sugar moiety) found"

    # Check for ester linkage between sugar hydroxyl and phosphate group
    # Look for an oxygen atom connected to both a sugar carbon and phosphorus
    ester_found = False
    for phosphate_match in phosphate_matches:
        # Get phosphorus atom index
        p_idx = phosphate_match[1]  # Phosphorus is the second atom in the pattern
        p_atom = mol.GetAtomWithIdx(p_idx)
        # Get oxygen atoms bonded to phosphorus
        phosphate_oxygen_idxs = [nbr.GetIdx() for nbr in p_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
        for sugar_match in sugar_matches:
            # Get atoms in sugar ring
            sugar_atom_idxs = list(sugar_match)
            # Check for ester linkage
            for idx in sugar_atom_idxs:
                atom = mol.GetAtomWithIdx(idx)
                # Look for oxygen atoms connected to this carbon
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 8:
                        nbr_idx = nbr.GetIdx()
                        # Check if this oxygen is connected to phosphorus
                        if nbr_idx in phosphate_oxygen_idxs:
                            ester_found = True
                            return True, "Contains monosaccharide with an alcoholic hydroxy group esterified with phosphoric acid"

    if not ester_found:
        return False, "Phosphate group not esterified to monosaccharide hydroxyl group"