"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
from rdkit import Chem

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is characterized by the presence of hydroperoxy groups (-OOH) 
    attached to a lipid-like structure, typically devoid of structural complexities like 
    charged groups or extensive cyclic elements.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for hydroperoxy group pattern (-OOH)
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2H]O")
    hydroperoxy_matches = mol.GetSubstructMatches(hydroperoxy_pattern)
    if not hydroperoxy_matches:
        return False, "No hydroperoxy groups found"

    # Consider minimum reasonable carbon atoms threshold to define a lipid
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 16:  # Increased to 16 for better alignment with lipid length
        return False, f"Only {c_count} carbon atoms found, potentially too few for lipid structure"

    # Check for presence of carboxylic group and avoid negative charges
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[O]")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_pattern)
    if not carboxylic_matches:
        return False, "No carboxylic structures found"
    
    # Exclude charged species
    if any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms()):
        return False, "Charged species present, excluding for lipid hydroperoxide"

    # Avoid extensive cyclic structures or complex multi-oxygen linkages (heuristically)
    # Such features commonly denote deviation from lipid hydroperoxides
    if mol.GetRingInfo().NumRings() > 1:
        return False, "Multiple rings detected, uncharacteristic of typical lipid hydroperoxides"

    return True, f"Contains {len(hydroperoxy_matches)} hydroperoxy group(s) and lipid-like structure with {c_count} carbon atoms"