"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
"""
Classifies: CDP-diacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    A CDP-diacylglycerol is a glycerol backbone with two fatty acid chains attached via ester bonds
    at positions 1 and 2, and a CDP (cytidine diphosphate) group attached via a phosphodiester bond
    at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a CDP-diacylglycerol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone with ester groups at positions 1 and 2
    # Glycerol backbone with two esterified hydroxyls at positions 1 and 2
    glycerol_pattern = Chem.MolFromSmarts("C(COC(=O)*)COC(=O)*")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No diacylglycerol backbone found"

    # Check for phosphodiester bond at position 3
    phosphodiester_pattern = Chem.MolFromSmarts("COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H](N2C=CC(=NC2=O)N)[C@@H](O)[C@H]1O")  # Simplified pattern for CDP group
    phosphodiester_matches = mol.GetSubstructMatches(phosphodiester_pattern)
    if not phosphodiester_matches:
        return False, "No CDP group attached via phosphodiester bond found"

    # Check for cytidine group attached to phosphate
    cytidine_pattern = Chem.MolFromSmarts("N2C=CC(=NC2=O)N")  # Pattern for cytosine base
    cytidine_matches = mol.GetSubstructMatches(cytidine_pattern)
    if not cytidine_matches:
        return False, "No cytidine (cytosine base) group found"

    # Check for fatty acid chains (long carbon chains attached via ester bonds)
    # Fatty acid chain pattern: Long aliphatic chain connected via ester linkage
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)OCC(C)C")  # Simplified pattern
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Found {len(fatty_acid_matches)} fatty acid chains, need at least 2"

    # Check molecular weight - CDP-diacylglycerol typically > 700 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, "Molecular weight too low for CDP-diacylglycerol"

    return True, "Contains glycerol backbone with two fatty acid chains and CDP group attached via phosphodiester bond"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17963',
        'name': 'CDP-diacylglycerol',
        'definition': 'A CDP-glycerol having unspecified acyl groups (most commonly fatty acyl groups) at the 1- and 2-positions.',
        'parents': []
    },
}