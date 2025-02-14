"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
"""
Classifies: CHEBI:87871 3'-hydroxyflavanone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    A 3'-hydroxyflavanone is a flavanone with a hydroxy substituent at position 3' of the phenyl ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for flavanone scaffold (C=C-C-C=O-C=C)
    flavanone_pattern = Chem.MolFromSmarts("[C&R1]=1[C&R1]=C[C&R1]([C&R1]=C[C&R1]=1)[C&R2]=2[C&R1]=C[C&R1]=C([C&R1]=2)[O&R1]")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone scaffold found"
    
    # Look for hydroxy group at 3' position of phenyl ring
    hydroxy_pattern = Chem.MolFromSmarts("[c&R1]1[c&R1][c&R1]([O&R1])[c&R1][c&R1][c&R1]1")
    hydroxy_match = mol.GetSubstructMatches(hydroxy_pattern)
    if not hydroxy_match:
        return False, "No hydroxy group at 3' position"
    
    # Check if the matched hydroxy group is part of the phenyl ring
    phenyl_ring_atoms = set(AllChem.GetMolFragFromAtomSmilesPattern(mol, "c1ccccc1", atomDepictionOptions=Chem.AtomDepictionOptions()))
    hydroxy_atom = hydroxy_match[0][2]
    if hydroxy_atom not in phenyl_ring_atoms:
        return False, "Hydroxy group not part of phenyl ring"
    
    return True, "Contains flavanone scaffold with hydroxy group at 3' position of phenyl ring"