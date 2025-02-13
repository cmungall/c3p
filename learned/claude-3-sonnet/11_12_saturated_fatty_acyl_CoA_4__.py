"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 11,12-saturated fatty acyl-CoA(4-) based on its SMILES string.
    A 11,12-saturated fatty acyl-CoA(4-) is a fatty acyl-CoA(4-) where the 11-12 bond of the fatty acyl group is saturated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 11,12-saturated fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for CoA substructure
    coa_pattern = Chem.MolFromSmarts("C(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)([O-])[O-])n1cnc2c(N)ncnc12)C(=O)NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA substructure found"
    
    # Check for fatty acyl chain
    acyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if not acyl_matches:
        return False, "No fatty acyl chain found"
    
    # Check for saturation at 11-12 position
    for acyl_match in acyl_matches:
        acyl_atoms = [mol.GetAtomWithIdx(idx) for idx in acyl_match]
        if len(acyl_atoms) < 12:
            continue  # Acyl chain too short
        
        atom11 = acyl_atoms[10]
        atom12 = acyl_atoms[11]
        if atom11.GetIsAromatic() or atom12.GetIsAromatic():
            continue  # Aromatic atoms
        
        bond = mol.GetBondBetweenAtoms(atom11.GetIdx(), atom12.GetIdx())
        if bond.GetBondType() == Chem.BondType.SINGLE:
            return True, "Fatty acyl chain has a saturated bond at position 11-12"
    
    return False, "Fatty acyl chain does not have a saturated bond at position 11-12"