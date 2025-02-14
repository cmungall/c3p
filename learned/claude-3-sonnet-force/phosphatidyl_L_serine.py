"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
Here's a program to classify if a given SMILES string represents a phosphatidyl-L-serine molecule:

"""
Classifies: CHEBI:17839 phosphatidyl-L-serine

A class of aminophospholipids in which a phosphatidyl group is esterified to the hydroxy group of serine.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a phosphatidyl-L-serine based on its SMILES string.

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
    
    # Look for phosphatidyl group (P-O-C-C-O-C(=O)-)
    phosphatidyl_pattern = Chem.MolFromSmarts("P(OCC(OC(=O))OC(=O))")
    if not mol.HasSubstructMatch(phosphatidyl_pattern):
        return False, "No phosphatidyl group found"
    
    # Look for serine residue (-CH(NH2)-C(=O)-O-)
    serine_pattern = Chem.MolFromSmarts("C(N)C(O)=O")
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "No serine residue found"
    
    # Look for ester linkage between phosphatidyl and serine
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, "Missing ester linkage between phosphatidyl and serine"
    
    # Look for fatty acid chains (long carbon chains)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, "Missing fatty acid chains"
    
    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"
    
    # Check molecular weight - phosphatidyl-L-serines typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for phosphatidyl-L-serine"
    
    return True, "Contains phosphatidyl group esterified to serine residue with fatty acid chains"