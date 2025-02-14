"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: CHEBI:36939 monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    A monoterpenoid is any terpenoid derived from a monoterpene, including compounds
    where the C10 skeleton has been rearranged or modified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for isoprenoid units
    ipr_pattern = Chem.MolFromSmarts("[C@H]1CC=C(C)C1")
    has_ipr_unit = mol.HasSubstructMatch(ipr_pattern)

    # Check for common monoterpene skeletons
    pinane_pattern = Chem.MolFromSmarts("[C@@]12C[C@H]([C@@H]([C@H]1[C@@H]3C[C@H](C3)C2)C)C")
    fenchane_pattern = Chem.MolFromSmarts("[C@H]12C[C@@H](C[C@H]1C2)C")
    has_mono_skeleton = mol.HasSubstructMatch(pinane_pattern) or mol.HasSubstructMatch(fenchane_pattern)

    # Check for common monoterpene functional groups
    alcohol_pattern = Chem.MolFromSmarts("OC")
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ether_pattern = Chem.MolFromSmarts("CO")
    aldehyde_pattern = Chem.MolFromSmarts("C=O")
    ketone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX3]")
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    sulfur_pattern = Chem.MolFromSmarts("[#16]")
    has_mono_groups = any([mol.HasSubstructMatch(p) for p in [alcohol_pattern, ester_pattern, ether_pattern,
                                                            aldehyde_pattern, ketone_pattern, acid_pattern,
                                                            sulfur_pattern]])

    # Check molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    if c_count == 10 and mol_wt < 300 and n_rings <= 2 and n_rotatable <= 5:
        if has_ipr_unit or has_mono_skeleton or has_mono_groups:
            return True, "Molecule contains features characteristic of monoterpenoids"
        else:
            return False, "Molecule lacks typical monoterpenoid features"
    else:
        return False, "Molecule does not meet basic criteria for monoterpenoids"