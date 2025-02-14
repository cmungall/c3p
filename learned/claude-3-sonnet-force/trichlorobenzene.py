"""
Classifies: CHEBI:27096 trichlorobenzene
"""
"""
Classifies: CHEBI:32221 trichlorobenzene

Any member of the class of chlorobenzenes carrying three chloro substituents at unspecified positions.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trichlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count number of chlorine atoms
    chlorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17)
    if chlorine_count != 3:
        return False, f"Molecule contains {chlorine_count} chlorine atoms, should be 3"
    
    # Check if there is a benzene ring
    benzene_ring = mol.GetRingInfo().AtomRings()
    has_benzene = any(len(ring) == 6 and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring) for ring in benzene_ring)
    if not has_benzene:
        return False, "No benzene ring found"
    
    # Check if all chlorines are attached to the benzene ring
    benzene_atoms = set()
    for ring in benzene_ring:
        if len(ring) == 6:
            benzene_atoms.update(ring)
            break
    
    chlorine_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17]
    if not all(any(mol.GetBondBetweenAtoms(i, j).GetIsAromatic() for j in benzene_atoms) for i in chlorine_atoms):
        return False, "Not all chlorine atoms are attached to the benzene ring"
    
    return True, "Molecule contains a benzene ring with three chlorine substituents"