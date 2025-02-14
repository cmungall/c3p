"""
Classifies: CHEBI:17334 penicillin
"""
"""
Classifies: penicillin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    A penicillin is defined as any substituted penam containing two methyl groups at position 2,
    a carboxylate group at position 3, and a carboxamido group at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define beta-lactam ring pattern (4-membered cyclic amide)
    beta_lactam = Chem.MolFromSmarts('C1CNC1=O')
    if not mol.HasSubstructMatch(beta_lactam):
        return False, "Beta-lactam ring not found"

    # Define thiazolidine ring pattern (5-membered ring with S and N)
    thiazolidine = Chem.MolFromSmarts('C1CSCN1')
    if not mol.HasSubstructMatch(thiazolidine):
        return False, "Thiazolidine ring not found"

    # Check if beta-lactam and thiazolidine rings are fused
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    fused = False
    for ring1 in atom_rings:
        if len(ring1) == 4:  # beta-lactam ring
            for ring2 in atom_rings:
                if len(ring2) == 5:  # thiazolidine ring
                    shared_atoms = set(ring1) & set(ring2)
                    if len(shared_atoms) >= 1:
                        fused = True
                        break
            if fused:
                break
    if not fused:
        return False, "Beta-lactam and thiazolidine rings are not fused"

    # Check for two methyl groups at position 2 (carbon adjacent to sulfur)
    methylated_carbon = Chem.MolFromSmarts('[S][C](C)(C)')
    if not mol.HasSubstructMatch(methylated_carbon):
        return False, "Two methyl substituents at position 2 not found"

    # Check for carboxylate group at position 3
    carboxylate = Chem.MolFromSmarts('C(=O)[O-,O][C]')
    if not mol.HasSubstructMatch(carboxylate):
        return False, "Carboxylate group at position 3 not found"

    # Check for carboxamido group at position 6
    carboxamido = Chem.MolFromSmarts('NC(=O)[C]')
    if not mol.HasSubstructMatch(carboxamido):
        return False, "Carboxamido group at position 6 not found"

    return True, "Molecule matches penicillin structure with required substituents"

__metadata__ = {
    'chemical_class': {
        'name': 'penicillin',
        'definition': 'Any member of the group of substituted penams containing two methyl substituents at position 2, a carboxylate substituent at position 3 and a carboxamido group at position 6.'
    }
}