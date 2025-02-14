"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: CHEBI:64381 3-sn-phosphatidyl-L-serine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    A 3-sn-phosphatidyl-L-serine has a glycerol backbone with acyl groups at positions 1 and 2,
    and a phospho-L-serine group at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-sn-phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define glycerol backbone with correct stereochemistry (sn-3 configuration)
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](CO[P](=O)(O)OC[C@H](N)C(=O)O)O")  # sn-3 glycerol backbone
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No sn-3 glycerol backbone with phospho-L-serine found"
    
    # Define ester bonds at positions 1 and 2
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")  # Ester group
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"
    
    # Check for phospho-L-serine group
    phospho_serine_pattern = Chem.MolFromSmarts("O[P](=O)(O)OC[C@H](N)C(=O)O")  # Phospho-L-serine
    if not mol.HasSubstructMatch(phospho_serine_pattern):
        return False, "Phospho-L-serine group not found"
    
    # Check for fatty acid chains (long alkyl chains attached via esters)
    fatty_acid_pattern = Chem.MolFromSmarts("OC(=O)C([CH2])CCCC")  # Simplified pattern for fatty acids
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Found {len(fatty_acid_matches)} fatty acid chains, need at least 2"
    
    # Check chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    required_centers = [atom.GetIdx() for atom in mol.GetAtoms() if atom.HasProp('_ChiralityPossible')]
    if len(chiral_centers) < len(required_centers):
        return False, "Chiral centers do not match required stereochemistry"
    
    return True, "Molecule is a 3-sn-phosphatidyl-L-serine with correct glycerol backbone and substituents"