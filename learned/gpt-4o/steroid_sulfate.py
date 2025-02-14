"""
Classifies: CHEBI:16158 steroid sulfate
"""
from rdkit import Chem

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    A steroid sulfate is a sulfuric ester of a steroid, featuring a steroid backbone
    with an attached sulfate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid sulfate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for sulfuric ester group (e.g., -OS(=O)(=O)=O)
    sulfate_group_pattern = Chem.MolFromSmarts("OS(=O)(=O)=O")
    if not mol.HasSubstructMatch(sulfate_group_pattern):
        return False, "No sulfate group found"
    
    # Identify steroid backbone pattern (sterane core structure C15H24)
    steroid_core_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4CC(C2)(C3)CC4C1")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core structure found"

    # Verify sulfate group's oxygen is connected to a steroid-like structure
    sulfate_attachments = mol.GetSubstructMatches(sulfate_group_pattern)
    attachment_found = False
    for match in sulfate_attachments:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Check if connected to any sterane carbon atoms
            for neighbor in atom.GetNeighbors():
                if neighbor.IsInRing() and neighbor.GetAtomicNum() == 6:
                    attachment_found = True
                    break
        if attachment_found:
            break

    if not attachment_found:
        return False, "Sulfate group not attached to steroid structure"
    
    return True, "Contains steroid structure with sulfate group attached via oxygen"

# Example usage (test with one of the SMILES from the list above):
print(is_steroid_sulfate("C1=C(OS(O)(=O)=O)C=CC2=C1CC[C@]3([C@@]4(CC[C@]([C@]4(CC[C@@]32[H])C)(C#C)O)[H])[H]"))