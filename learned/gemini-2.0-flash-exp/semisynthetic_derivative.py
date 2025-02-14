"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on its SMILES string.
    This version focuses on identifying a complex scaffold and synthetic modifications.
    It uses SMARTS patterns for common modifications and attempts a core identification

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a semisynthetic derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    num_heavy_atoms = mol.GetNumHeavyAtoms()
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    num_chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    
    if num_heavy_atoms < 10: #Added min heavy atoms
        return False, "Too few heavy atoms to be a semisynthetic derivative."
    
    if num_rings <= 1 and num_chiral_centers < 2:
        #Initial check to filter linear/simple structures
         return False, "Too simple, lacks ring structures and chiral centers"

    # Attempt to identify a complex core
    # Looking for a complex core with rings
    
    # Pattern for 3 or more rings connected by single bonds
    complex_ring_pattern1 = Chem.MolFromSmarts("[R]~[R]~[R]")
    # Pattern for 2 or more rings fused to a non-aromatic ring.
    complex_ring_pattern2 = Chem.MolFromSmarts("[R;!a]1@[R]@[R]1")
    
    has_complex_core = mol.HasSubstructMatch(complex_ring_pattern1) or mol.HasSubstructMatch(complex_ring_pattern2)

    if not has_complex_core:
        return False, "Molecule lacks a complex core, unlikely a semi-synthetic."
   

    # Common synthetic modifications (esters, amides, ethers, halogens)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    num_esters = len(mol.GetSubstructMatches(ester_pattern))
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    num_amides = len(mol.GetSubstructMatches(amide_pattern))
    ether_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    num_ethers = len(mol.GetSubstructMatches(ether_pattern))
    num_halogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [9, 17, 35, 53])
    
    modification_score = 0
    modification_score += num_esters * 1 # Ester are common
    modification_score += num_amides * 2 # Amides are less common
    modification_score += num_ethers * 1  # Ethers common
    
    # Penalize excessive simple modifications
    if (num_esters + num_amides + num_ethers) > 5:
        modification_score = modification_score -2

    if num_halogens > 0: # Halogens should increase score, but not too much
         modification_score = modification_score + 1

    #Combined check
    if modification_score > 0:
            return True, "Contains structural features suggestive of a semisynthetic derivative"
    else:
        return False, "No significant synthetic modifications found"