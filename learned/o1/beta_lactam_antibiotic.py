"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: CHEBI:35627 beta-lactam antibiotic
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    A beta-lactam antibiotic is an organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define beta-lactam ring pattern
    beta_lactam_smarts = '[N;R][C;R][C;R][C;R](=O)'
    beta_lactam_pattern = Chem.MolFromSmarts(beta_lactam_smarts)
    
    # Check for beta-lactam ring
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"
    
    # Check for organonitrogen heterocyclic nature
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    is_organonitrogen_heterocycle = False
    
    # Iterate over rings to find rings containing nitrogen
    for ring in atom_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Check if any atom in the ring is nitrogen
        if any(atom.GetAtomicNum() == 7 for atom in ring_atoms):
            is_organonitrogen_heterocycle = True
            break
    
    if not is_organonitrogen_heterocycle:
        return False, "Molecule is not an organonitrogen heterocycle"
    
    return True, "Contains beta-lactam ring and is an organonitrogen heterocycle"