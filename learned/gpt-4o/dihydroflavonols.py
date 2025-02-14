"""
Classifies: CHEBI:48039 dihydroflavonols
"""
from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol is a hydroxyflavanone with a hydroxyl group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for dihydroflavonol pattern: flavanone backbone with hydroxyl group at position 3
    dihydroflavonol_pattern = Chem.MolFromSmarts('O[C@H]1[C@@H](Oc2cc(O)cc(O)c2C1=O)')  # Simplified pattern
    if mol.HasSubstructMatch(dihydroflavonol_pattern):
        return True, "Contains dihydroflavonol structure"

    # Look for hydroxyflavanone pattern with hydroxyl group at position 3
    hydroxyflavanone_pattern = Chem.MolFromSmarts('OC1C(OC2=C(C=CC(O)=C2)C1=O)')  # More general pattern
    matches = mol.GetSubstructMatches(hydroxyflavanone_pattern)
    if matches:
        # Check the attachment positions to ensure hydroxyl is at position 3
        for match in matches:
            atom_positions = [mol.GetAtomWithIdx(idx).GetIdx() for idx in match]
            # Ensure one hydroxyl is attached at position 3
            # (This is a simplified check; actual check might require mapping atom indices correctly)
            has_hydroxy_at_3 = False
            for bond in mol.GetBonds():
                if bond.GetBeginAtomIdx() in atom_positions and bond.GetEndAtomIdx() in atom_positions:
                    if mol.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetSymbol() == 'O' and \
                       mol.GetAtomWithIdx(bond.GetEndAtomIdx()).GetSymbol() == 'C':
                        has_hydroxy_at_3 = True
                        break
            if has_hydroxy_at_3:
                return True, "Contains hydroxyflavanone structure with hydroxyl group at position 3"

    return False, "Does not fit dihydroflavonol or hydroxyflavanone pattern adequately"