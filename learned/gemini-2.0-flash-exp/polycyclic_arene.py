"""
Classifies: CHEBI:33848 polycyclic arene
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene based on its SMILES string.
    A polycyclic arene consists of fused aromatic rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polycyclic arene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of aromatic rings
    aromatic_ring_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(aromatic_ring_pattern):
        return False, "No aromatic rings found"
    
    # Check if there are multiple aromatic rings
    aromatic_count = len(mol.GetSubstructMatches(aromatic_ring_pattern))
    if aromatic_count < 2:
         return False, "Less than 2 aromatic rings found, not a polycyclic arene"
    
    #check presence of fused rings
    fused_ring_pattern = Chem.MolFromSmarts("c12ccccc1cccc2")
    if not mol.HasSubstructMatch(fused_ring_pattern):
        return False, "No fused rings detected"
    
    # Check for non-aromatic atoms other than hydrogens; this part is to eliminate complex molecules with other functional groups.
    non_aromatic_atoms = [atom for atom in mol.GetAtoms() if atom.GetIsAromatic() == False and atom.GetAtomicNum() != 1]
    if len(non_aromatic_atoms) > 0:
        return False, f"Non-aromatic atoms other than H found {len(non_aromatic_atoms)}"

    # Check for single bonds between aromatic rings to exclude cases where they are just joined not fused
    single_bond_between_aromatic_rings = Chem.MolFromSmarts("[c]~[c]")
    if len(mol.GetSubstructMatches(single_bond_between_aromatic_rings)) == 0:
            return False, "Rings are not fused"
        

    # Check that at least 50% of the atoms are aromatic.
    total_atoms = mol.GetNumAtoms()
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())

    if total_atoms == 0:
        return False, "Molecule has no atoms"

    if aromatic_atoms/total_atoms < 0.5:
        return False, "Not enough aromatic atoms"



    return True, "Contains fused aromatic rings, classified as polycyclic arene"