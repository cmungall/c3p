"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: clavulone
Definition: A class of esterified prostanoids obtained from marine corals.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a clavulone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core cyclopentenone patterns
    # More flexible patterns to catch different substitution variants
    core_patterns = [
        # Basic cyclopentenone core with various substitutions
        Chem.MolFromSmarts("[#6]1[#6](=[O])[#6]=[#6][#6]1"),
        # Halogenated variant
        Chem.MolFromSmarts("[Cl,Br,I][#6]1[#6](=[O])[#6][#6][#6]1"),
        # Conjugated variant
        Chem.MolFromSmarts("[#6]1[#6](=[O])[#6](=[#6])[#6][#6]1")
    ]
    
    # Check for at least one core pattern
    core_found = False
    for pattern in core_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            core_found = True
            break
    
    if not core_found:
        return False, "Missing cyclopentenone core structure"

    # Side chain patterns
    side_chain_patterns = [
        # Unsaturated aliphatic chain
        Chem.MolFromSmarts("C/C=C/CC"),
        # Ester-containing chain
        Chem.MolFromSmarts("CC(=O)OC"),
        # Conjugated unsaturated chain
        Chem.MolFromSmarts("C=CC=C"),
        # Acetoxy group
        Chem.MolFromSmarts("OC(=O)C")
    ]
    
    # Count matching side chain features
    side_chain_count = 0
    for pattern in side_chain_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            side_chain_count += 1
    
    if side_chain_count < 2:
        return False, "Insufficient characteristic side chain features"

    # Count key atoms and features
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    hal_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [9,17,35,53])
    ring_count = len(Chem.GetSymmSSSR(mol))
    
    # Basic requirements
    if not (15 <= c_count <= 35):
        return False, "Carbon count outside typical range for clavulones"
    
    if o_count < 2:
        return False, "Insufficient oxygen content"
        
    if ring_count != 1:
        return False, "Must contain exactly one ring"
        
    # Check for required functional groups
    ester_pattern = Chem.MolFromSmarts("C(=O)O[#6]")
    if ester_pattern is not None:
        ester_count = len(mol.GetSubstructMatches(ester_pattern))
        if ester_count < 1:
            return False, "Missing required ester group"

    # Additional check for conjugated system
    conjugated_pattern = Chem.MolFromSmarts("C=CC=C")
    if conjugated_pattern is not None and not mol.HasSubstructMatch(conjugated_pattern):
        return False, "Missing conjugated system"

    # Must have either halogen or multiple oxygen-containing groups
    if hal_count == 0 and o_count < 4:
        return False, "Insufficient characteristic substituents"

    return True, "Contains characteristic clavulone structural features including cyclopentenone core and required substitution patterns"