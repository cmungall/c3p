"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: hydroxynaphthoquinone
Definition: Any naphthoquinone in which the naphthaoquinone moiety is substituted by at least one hydroxy group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a hydroxynaphthoquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First check for naphthoquinone core
    # Two patterns to catch both ortho and para quinone arrangements
    naphthoquinone_pattern1 = Chem.MolFromSmarts("[#6]1[#6][#6][#6]2[#6](=[O])[#6][#6](=[O])[#6][#6]12") # ortho
    naphthoquinone_pattern2 = Chem.MolFromSmarts("[#6]1[#6][#6][#6]2[#6](=[O])[#6][#6][#6](=[O])[#6]12") # para
    
    has_naphthoquinone = mol.HasSubstructMatch(naphthoquinone_pattern1) or \
                         mol.HasSubstructMatch(naphthoquinone_pattern2)
    
    if not has_naphthoquinone:
        return False, "No naphthoquinone core found"
    
    # Get the atoms that are part of the naphthoquinone core
    core_atoms = set()
    if mol.HasSubstructMatch(naphthoquinone_pattern1):
        core_atoms.update(mol.GetSubstructMatch(naphthoquinone_pattern1))
    if mol.HasSubstructMatch(naphthoquinone_pattern2):
        core_atoms.update(mol.GetSubstructMatch(naphthoquinone_pattern2))
    
    # Look for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    if not hydroxyl_matches:
        return False, "No hydroxyl groups found"
    
    # Check if any hydroxyl group is attached to the naphthoquinone core
    for match in hydroxyl_matches:
        oh_oxygen = mol.GetAtomWithIdx(match[0])
        for neighbor in oh_oxygen.GetNeighbors():
            if neighbor.GetIdx() in core_atoms:
                position = "ortho" if mol.HasSubstructMatch(naphthoquinone_pattern1) else "para"
                return True, f"Contains {position}-naphthoquinone core with at least one hydroxyl group attached"
    
    return False, "No hydroxyl groups attached to naphthoquinone core"

def test_examples():
    """Test function with some example molecules"""
    examples = [
        "Oc1cc(O)c2C(=O)C=C(O)C(=O)c2c1",  # flaviolin
        "OC1=C([C@H]2CC[C@@H](CC2)c2ccc(Cl)cc2)C(=O)c2ccccc2C1=O",  # atovaquone
        "Oc1cccc2C(=O)C=CC(=O)c12",  # juglone
        "CC1=CC(=O)C=CC1=O"  # not a hydroxynaphthoquinone (methylquinone)
    ]
    
    for smi in examples:
        result, reason = is_hydroxynaphthoquinone(smi)
        print(f"SMILES: {smi}")
        print(f"Is hydroxynaphthoquinone: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_examples()