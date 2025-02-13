"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:36349 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.
    A long-chain fatty acyl-CoA(4-) is a fatty acyl-CoA with a long fatty acid chain (typically >14 carbons)
    and deprotonated phosphate and diphosphate groups (carrying a -4 charge).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA backbone
    coa_pattern = Chem.MolFromSmarts("C(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)([O-])[O-])n1cnc2c(N)ncnc12)C(=O)NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA backbone"
    
    # Look for fatty acyl chain (>14 carbons)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "Fatty acid chain too short (<14 carbons)"
    
    # Check for ester linkage between fatty acid and CoA
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, "Missing or multiple ester linkages"
    
    # Count negative charges (should be -4)
    charged_atoms = [atom for atom in mol.GetAtoms() if atom.GetFormalCharge() != 0]
    total_charge = sum(atom.GetFormalCharge() for atom in charged_atoms)
    if total_charge != -4:
        return False, f"Total charge is {total_charge}, should be -4"
    
    return True, "Contains a long fatty acid chain linked to the CoA backbone via an ester bond, with deprotonated phosphate and diphosphate groups"