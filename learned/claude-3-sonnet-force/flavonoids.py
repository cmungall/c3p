"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: CHEBI:27558 flavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid is a phenyl-substituted 1-phenylpropane derivative with a C15 or C16 skeleton,
    or a condensed structure with a C6-C3 lignan precursor. It may include flavonoids, isoflavonoids,
    neoflavonoids, chalcones, dihydrochalcones, aurones, pterocarpans, coumestarins, rotenoids,
    flavonolignans, homoflavonoids, and flavonoid oligomers.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for flavonoid skeleton patterns (C6-C3-C6, C6-C3-C3, or C6-C3-C6-C3)
    flavonoid_pattern = Chem.MolFromSmarts("c1ccc(cc1)cccc2ccc(cc2)Oc3ccccc3")  # C6-C3-C6
    isoflavonoid_pattern = Chem.MolFromSmarts("c1ccc(cc1)C=2C(=O)Oc3ccccc3C2=O")  # C6-C3-C3
    neoflavonoid_pattern = Chem.MolFromSmarts("c1ccc(cc1)cccc2cccc3ccccc32")  # C6-C3-C6-C3
    if not (mol.HasSubstructMatch(flavonoid_pattern) or
            mol.HasSubstructMatch(isoflavonoid_pattern) or
            mol.HasSubstructMatch(neoflavonoid_pattern)):
        return False, "No flavonoid skeleton found"

    # Count aromatic rings and check for condensed systems
    ring_info = mol.GetRingInfo()
    aromatic_rings = [r for r in ring_info.AtomRings() if ring_info.IsCycleAromaticHs(r)]
    n_aromatic_rings = len(aromatic_rings)
    if n_aromatic_rings > 3:
        return True, "Contains condensed flavonoid skeleton with more than 3 aromatic rings"

    # Check molecular weight - flavonoids typically 200-800 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 800:
        return False, "Molecular weight outside typical range for flavonoids"

    # Count oxygen atoms (common in flavonoids)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, "Too few oxygen atoms for flavonoid"

    return True, "Contains flavonoid skeleton with phenyl-substituted 1-phenylpropane derivative"