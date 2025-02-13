"""
Classifies: CHEBI:75874 para-terphenyl
"""
"""
Classifies: CHEBI:62648 para-terphenyl
A ring assembly based on a 1,4-diphenylbenzene skeleton and its substituted derivatives thereof.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a para-terphenyl, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for 1,4-diphenylbenzene core
    core_pattern = Chem.MolFromSmarts("c1ccc(cc1)-c2ccc(cc2)-c3ccccc3")
    match = mol.GetSubstructMatches(core_pattern)
    if not match:
        return False, "Does not contain 1,4-diphenylbenzene core"
    
    # Check for aromatic rings and substituents
    rings = mol.GetRingInfo().AtomRings()
    aromatic_rings = [r for r in rings if all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in r)]
    if len(aromatic_rings) < 3:
        return False, "Does not have at least three aromatic rings"
    
    substituted_rings = []
    for ring in aromatic_rings:
        ring_atoms = [mol.GetAtomWithIdx(a) for a in ring]
        if any(a.GetTotalNumHs() < 1 for a in ring_atoms):
            substituted_rings.append(ring)
    if len(substituted_rings) < 3:
        return False, "Core rings are not substituted"
    
    return True, "Contains 1,4-diphenylbenzene core with substituted derivatives"