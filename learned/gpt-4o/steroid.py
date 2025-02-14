"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid is identified by the presence of a cyclopenta[a]phenanthrene carbon skeleton,
    typically comprising three 6-membered rings and one 5-membered ring, with flexibility for common
    functional modifications.

    Args:
        smiles (str): SMILES string of the chemical compound

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Examine ring structures - steroids have three 6-membered rings and one 5-membered ring
    ring_info = mol.GetRingInfo()
    sssr = ring_info.NumRings()
    
    if sssr < 4:
        return False, "Fewer than 4 rings identified, not a typical steroid ring system"
    
    six_membered_rings = sum(1 for ring in ring_info.AtomRings() if len(ring) == 6)
    five_membered_rings = sum(1 for ring in ring_info.AtomRings() if len(ring) == 5)
    
    if six_membered_rings < 3 or five_membered_rings < 1:
        return False, "Steroid structure not found in terms of 6:5 ring pattern"
    
    # Check for presence of typical functional groups and common modifications
    # Example: methyl groups at specific positions which are common in steroids
    methyl_smarts = Chem.MolFromSmarts("[C](C)(C)[CH](C)")
    if not mol.HasSubstructMatch(methyl_smarts):
        return False, "Missing typical methyl groups seen in many steroids"
    
    # Consider typical molecular weight range of steroids if needed
    mol_weight = Descriptors.MolWt(mol)
    if mol_weight < 250 or mol_weight > 500:
        return False, "Molecular weight atypical for steroids"
    
    # Optionally, add more checks for specific functional groups (hydroxyl, carbonyl, etc.)
    # Hydroxyl presence (common in many steroids)
    hydroxyl_smarts = Chem.MolFromSmarts("[OX2H]")
    if mol.HasSubstructMatch(hydroxyl_smarts):
        return True, "Molecule matches steroid structure with hydroxyl group."
    
    return True, "Molecule matches typical steroid structure"

# Examples SMILES test (these should return True)
smiles_examples = [
    "C[C@]12CC[C@H]3[C@H]([C@H]1CC[C@]2(C#C)O)CCC4=C3C=CC(=C4)OC",  # Expected True
    "CC(C)=CCC[C@](C)(O)[C@H]1CC[C@]2(C)[C@@H]1[C@H](O)C[C@@H]1[C@@]3(C)CC[C@H](O)C(C)(C)[C@@H]3CC[C@@]21C",  # True
]

for example in smiles_examples:
    print(is_steroid(example))