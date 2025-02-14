"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    Aldoses are polyhydroxy aldehydes (H[CH(OH)]nC(=O)H, n >= 2) and their intramolecular hemiacetals.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldose, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon and oxygen atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)


    # Check for open-chain aldose structure
    # Look for a terminal aldehyde, followed by a chain of at least 2 carbons, all with OH groups
    open_chain_pattern = Chem.MolFromSmarts("[CH1](=O)[CH]([OX2H])[CH]([OX2H])([#6])[#6]")
    if mol.HasSubstructMatch(open_chain_pattern):
        
        # Verify the specific number of oxygens for open chain forms: n_carbons + 1
        if oxygen_count != carbon_count+1 :
            return False, "Open chain detected but not an aldose: wrong number of oxygens"
        return True, "Open-chain aldose structure detected"


    # Check for cyclic hemiacetal form (5 or 6 membered ring with an anomeric carbon)
    # Pattern for 5-membered ring
    hemiacetal_pattern5 = Chem.MolFromSmarts("[OX2]1[CX4H]([OX2H])[CX4H]([OX2H])[CX4H]([OX2H])[CX4H]1")
    # Pattern for 6-membered ring
    hemiacetal_pattern6 = Chem.MolFromSmarts("[OX2]1[CX4H]([OX2H])[CX4H]([OX2H])[CX4H]([OX2H])[CX4H]([OX2H])[CX4H]1")

    if mol.HasSubstructMatch(hemiacetal_pattern5) or mol.HasSubstructMatch(hemiacetal_pattern6):
         # Check the ring is formed only by carbon and oxygen
        ring_atoms_pattern5 = Chem.MolFromSmarts("[OX2]1[CX4H]([OX2H])[CX4H,OX2][CX4H,OX2][CX4H,OX2]1")
        ring_atoms_pattern6 = Chem.MolFromSmarts("[OX2]1[CX4H]([OX2H])[CX4H,OX2][CX4H,OX2][CX4H,OX2][CX4H,OX2]1")
        if not (mol.HasSubstructMatch(ring_atoms_pattern5) or mol.HasSubstructMatch(ring_atoms_pattern6)):
            return False, "Cyclic hemiacetal detected, but not a carbohydrate ring"
        
        # Check there are enough carbons
        if carbon_count < 3:
            return False, "Cyclic aldose must contain at least 3 carbon atoms"
        # Count the number of OH groups (excluding the anomeric one) in the ring
        ring_hydroxy_pattern = Chem.MolFromSmarts("[CX4H]([OX2H])")
        ring_matches = mol.GetSubstructMatches(ring_hydroxy_pattern)
        if len(ring_matches) < 2: #must have at least 2 hydroxyl groups, not counting the anomeric one
            return False, "Cyclic aldose must have at least 2 other OH groups"
        return True, "Cyclic hemiacetal aldose structure detected"


    return False, "Does not match aldose criteria"