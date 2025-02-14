"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    Diterpenoids have a C20 skeleton (or modified) and are built from isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Check for ring systems complexity. Use ring bond count
    ring_bond_count = rdMolDescriptors.CalcNumSpiroAtoms(mol) + rdMolDescriptors.CalcNumBridgeheadAtoms(mol) + rdMolDescriptors.CalcNumRings(mol)

    if ring_bond_count < 3:
        return False, f"Diterpenoids should have a complex ring system (>=3 ring bonds) this one has {ring_bond_count}"

    # Check for presence of multiple isoprene units - less strict pattern this time
    isoprene_pattern = Chem.MolFromSmarts("[CX4]([CX4])([CX4])~[CX4]~[CX4]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)

    if len(isoprene_matches) < 2:
        return False, "Too few isoprene units detected (less than 2)"
    
    if c_count < 18 or c_count > 22 :
         return False, f"Carbon count {c_count} is not within the diterpenoid range (18-22)"


    # Exclude common steroidal substructures - tetracyclic with specific substitutions and 17 carbons typically.
    steroid_pattern = Chem.MolFromSmarts("[C]1[C][C]2[C]3[C]([C]1)[C]([C]2)[C][C]([C]3)[C]")
    if mol.HasSubstructMatch(steroid_pattern) and c_count < 21 :
        return False, "Steroidal substructure detected, likely not a diterpenoid"

    return True, "Likely a diterpenoid based on ring system and multiple isoprene units"