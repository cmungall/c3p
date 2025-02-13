"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: CHEBI:28413 flavonols
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    Flavonols are hydroxyflavones with a hydroxy group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic flavonoid core with 3-OH group and 4-ketone (flavonol backbone)
    # [#6]1=[#6]-c2c([#6](=[O])[#6]([OH])=[#6]1)c([OH])cc([OH])c2
    flavonol_pattern = Chem.MolFromSmarts('[#6]1=[#6]-c2c([#6](=[O])[#6]([OH])=[#6]1)c([OH,O-1])cc([OH,O-1])c2')
    
    # Alternative pattern to catch variations
    flavonol_pattern2 = Chem.MolFromSmarts('O=C1C(O)=C(Oc2ccccc12)c1ccccc1')
    
    if not (mol.HasSubstructMatch(flavonol_pattern) or mol.HasSubstructMatch(flavonol_pattern2)):
        return False, "Missing flavonol core structure (3-hydroxyflavone backbone)"

    # Check for ketone at position 4
    ketone_pattern = Chem.MolFromSmarts('O=C1c2ccccc2OC(c2ccccc2)=C1O')
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "Missing ketone group at position 4"

    # Count oxygen atoms (should have at least 4: 3-OH, 4=O, and ring O)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, "Insufficient oxygen atoms for flavonol structure"

    # Check for aromatic rings (should have at least 2)
    aromatic_rings = 0
    for atom in mol.GetAtoms():
        if atom.IsInRing() and atom.GetIsAromatic():
            aromatic_rings += 1
    if aromatic_rings < 8:  # Each aromatic ring contributes multiple aromatic atoms
        return False, "Missing required aromatic rings"

    # Additional check for 2,3-double bond
    double_bond_pattern = Chem.MolFromSmarts('C1=COc2ccccc2C1=O')
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "Missing required double bond between positions 2 and 3"

    return True, "Contains 3-hydroxyflavone backbone with required structural features"