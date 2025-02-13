"""
Classifies: CHEBI:33913 corrinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_corrinoid(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    A corrinoid contains four reduced or partly reduced pyrrole rings joined in a 
    macrocycle by three =C- groups and one direct carbon-carbon bond linking alpha positions.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple[bool, str]: (is_corrinoid, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for corrin nucleus - a more specific pattern for the core structure
    # This pattern describes the four nitrogen atoms in their characteristic arrangement
    # with connecting carbons, allowing for different oxidation states and substituents
    corrin_core = Chem.MolFromSmarts('[#7]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]')
    if not mol.HasSubstructMatch(corrin_core):
        return False, "Missing corrin nucleus structure"

    # Look for the four modified pyrrole rings with more flexible matching
    # This pattern allows for reduced and substituted forms
    pyrrole_pattern = Chem.MolFromSmarts('[#7]([#6]~[#6])~[#6]~[#6]')
    matches = mol.GetSubstructMatches(pyrrole_pattern)
    if len(matches) < 4:
        return False, "Missing required pyrrole-like rings"

    # Check for presence of cobalt (common in corrinoids but not required)
    has_cobalt = any(atom.GetAtomicNum() == 27 for atom in mol.GetAtoms())
    
    # Count nitrogens in the core structure
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 4:
        return False, "Insufficient nitrogen atoms for corrin structure"

    # Check molecular size - corrinoids are large molecules
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 20:
        return False, "Molecule too small for corrinoid structure"

    # Look for characteristic carbon-carbon connections
    # This includes both double bonds and single bonds between rings
    ring_connections = Chem.MolFromSmarts('[#7]~[#6](~[#6])~[#6]~[#7]')
    if not mol.HasSubstructMatch(ring_connections):
        return False, "Missing required ring connections"

    # Success message varies depending on cobalt presence
    if has_cobalt:
        reason = ("Contains corrin nucleus with four modified pyrrole rings "
                 "and cobalt coordination center")
    else:
        reason = ("Contains corrin nucleus with four modified pyrrole rings "
                 "in characteristic arrangement")

    return True, reason