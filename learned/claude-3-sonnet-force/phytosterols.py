"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: CHEBI:37616 phytosterols
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Lipinski import NumRotatableBonds

def is_phytosterol(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols similar to cholesterol that occur in plants and vary only in carbon side chains and/or presence or absence of a double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the sterol backbone (tetracyclic triterpene)
    sterol_pattern = Chem.MolFromSmarts("[C@]1([C@H]([C@@H]2[C@@]([C@]([C@H](C1)CC2)(C)C)(C)C)C")
    if not mol.HasSubstructMatch(sterol_pattern):
        return False, "No sterol backbone found"
    
    # Check for hydroxyl group at position 3
    hydroxy_pattern = Chem.MolFromSmarts("[C@H](O)[C@H]1[C@@]2([C@H](C[C@@H]3[C@]([C@@](CC3)(C)[H])(C)[H])C2)")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No hydroxyl group at position 3"
    
    # Check for side chain length and flexibility
    n_rotatable = NumRotatableBonds(mol)
    if n_rotatable < 2:
        return False, "Side chain too rigid for phytosterol"
    
    # Check molecular weight - phytosterols typically ~400-500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 550:
        return False, "Molecular weight outside typical range for phytosterols"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 27 or c_count > 32:
        return False, "Incorrect carbon count for phytosterol"
    if o_count != 1:
        return False, "Must have exactly 1 oxygen (hydroxyl group)"
    
    return True, "Contains sterol backbone with hydroxyl at position 3 and flexible side chain"