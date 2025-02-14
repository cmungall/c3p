"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: CHEBI:33779 organic sulfide
Organic sulfides are compounds having the structure RSR (R =/= H)
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide based on its SMILES string.
    An organic sulfide has the structure RSR (R != H).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find sulfur atoms
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    
    # Check if there is at least one sulfur atom
    if not sulfur_atoms:
        return False, "No sulfur atoms found"
    
    # Exclude thiocyanate and thiol groups
    thiocyanate_pattern = Chem.MolFromSmarts("[S-]C#N")
    thiol_pattern = Chem.MolFromSmarts("[SH]")
    if mol.HasSubstructMatch(thiocyanate_pattern) or mol.HasSubstructMatch(thiol_pattern):
        return False, "Molecule contains a thiocyanate or thiol group"
    
    # Check if any sulfur atom satisfies the RSR structure
    for sulfur in sulfur_atoms:
        neighbors = [atom for bond in sulfur.GetBonds() for atom in [bond.GetBeginAtom(), bond.GetEndAtom()] if atom.GetAtomicNum() != 1]
        if len(neighbors) >= 2:
            r1, r2 = neighbors[:2]
            if r1.GetAtomicNum() != 1 and r2.GetAtomicNum() != 1:
                return True, "Contains the structure RSR where R is not hydrogen (organic sulfide)"
    
    return False, "No sulfur atom satisfies the RSR structure"