"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
"""
Classifies: phosphatidyl-L-serine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a phosphatidyl-L-serine based on its SMILES string.
    A phosphatidyl-L-serine is a glycerophospholipid with a phosphatidyl group esterified
    to the hydroxy group of L-serine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone (3 connected carbon atoms)
    glycerol_pattern = Chem.MolFromSmarts("[CH2][CH][CH2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for at least 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Check for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for serine moiety connected to phosphate
    # Serine pattern: N[C@H](C(=O)O)CO
    # Phosphate connected to serine hydroxy group
    serine_phosphate_pattern = Chem.MolFromSmarts("N[C@H](C(=O)O)CO[P](=O)(O)O")
    if not mol.HasSubstructMatch(serine_phosphate_pattern):
        return False, "No serine moiety connected to phosphate group found"

    # All criteria met
    return True, "Contains glycerol backbone with two fatty acid chains, phosphate group, and serine moiety"

__metadata__ = {
    'chemical_class': {
        'name': 'phosphatidyl-L-serine',
        'definition': 'A class of aminophospholipids in which a phosphatidyl group is esterified to the hydroxy group of serine.'
    }
}