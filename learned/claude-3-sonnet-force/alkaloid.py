"""
Classifies: CHEBI:22315 alkaloid
"""
"""
Classifies: CHEBI:3115 alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.
    Alkaloids are naturally occurring, basic nitrogen compounds (mostly heterocyclic),
    found in plants, bacteria, fungi, and animals. They typically contain a nitrogen
    atom in a heterocyclic ring system.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for nitrogen atom
    if not any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "No nitrogen atom found"

    # Check for heterocyclic rings
    ring_info = mol.GetRingInfo()
    hetero_rings = [ring for ring in ring_info.AtomRings() if any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in ring)]
    if not hetero_rings:
        return False, "No heterocyclic rings found"

    # Check for basic nitrogen (N-C bond)
    basic_n_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and any(mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx()).GetBondType() == Chem.BondType.SINGLE for nbr in atom.GetNeighbors())]
    if not basic_n_atoms:
        return False, "No basic nitrogen atoms found"

    # Check for common alkaloid heterocycles
    common_rings = ['pyridine', 'piperidine', 'quinoline', 'isoquinoline', 'indole', 'purine', 'imidazole', 'pyrazine', 'pyrimidine', 'pyrrole', 'pyrrolidine', 'tropane', 'pyrrolizidine']
    ring_smarts = [''.join(['n', str(len(ring)), 'aaannnnn']) for ring in hetero_rings]
    ring_patterns = [Chem.MolFromSmarts(smart) for smart in ring_smarts]
    ring_types = [AllChem.GetMolTemplateSetScore(mol, common_rings) for ring_pattern in ring_patterns]
    if any(ring_type > 0 for ring_type in ring_types):
        return True, "Contains nitrogen atom and heterocyclic ring(s) common in alkaloids"

    return True, "Contains nitrogen atom and heterocyclic ring(s)"