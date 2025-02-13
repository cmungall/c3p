"""
Classifies: CHEBI:47787 11-oxo steroid
"""
"""
Classifies: 11-oxo steroid
Definition: Any oxo steroid that has an oxo substituent at position 11
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (True if molecule is an 11-oxo steroid, reason for classification)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for basic steroid core (four fused rings)
    steroid_core = Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~[#6]~1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # SMARTS pattern for 11-oxo group in steroid context
    # This pattern looks for the specific environment of the 11-position ketone
    # in the steroid ring system
    oxo_11_pattern = Chem.MolFromSmarts("[#6]1~[#6]~[#6]2~[#6]~[#6]~[#6]~3~[#6](=[O:1])~[#6]~[#6]~[#6]~[#6]~3~[#6]~2~[#6]~[#6]~1")
    
    # Find matches for the 11-oxo pattern
    matches = mol.GetSubstructMatches(oxo_11_pattern)
    if not matches:
        return False, "No ketone group at position 11"
    
    # Count carbons to verify it's in the typical steroid range
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 19 or carbon_count > 30:
        return False, f"Carbon count ({carbon_count}) outside typical steroid range (19-30)"
        
    # Additional check for reasonable molecular weight
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 500:
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical steroid range (250-500)"
        
    # Count rings to ensure we have the right ring system
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"
        
    return True, "Molecule contains steroid core with ketone group at position 11"