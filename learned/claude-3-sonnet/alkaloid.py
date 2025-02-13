"""
Classifies: CHEBI:22315 alkaloid
"""
"""
Classifies: CHEBI:3230 alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.
    Alkaloids are naturally occurring, basic nitrogen compounds (mostly heterocyclic)
    found in plants, bacteria, fungi, and animals. They typically have one or more
    nitrogen-containing heterocyclic rings.

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
    
    # Check for presence of basic nitrogen
    has_basic_nitrogen = any(atom.GetAtomicNum() == 7 and atom.GetTotalNumHs() < 3 for atom in mol.GetAtoms())
    if not has_basic_nitrogen:
        return False, "No basic nitrogen present"
    
    # Check for heterocyclic rings
    rings = mol.GetRingInfo().AtomRings()
    has_heterocyclic_ring = any(any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in ring) for ring in rings)
    if not has_heterocyclic_ring:
        return False, "No heterocyclic rings found"
    
    # Check molecular weight - alkaloids typically 100-1000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 1000:
        return False, "Molecular weight outside typical range for alkaloids"
    
    # Count nitrogen atoms - alkaloids typically have 1-5
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 1 or n_count > 5:
        return False, "Number of nitrogen atoms outside typical range for alkaloids"
    
    # Check for common alkaloid substructures (optional)
    alkaloid_patterns = ['[NX3+]', '[NX4+]', '[N!X3]([!X1])!@[!X1]', '[N!X3]=,:[N]',
                         '[N!X3]=[N]', '[N!X3+](=[N])[!X1]', '[N!X3+]=C[N!X3]']
    has_alkaloid_pattern = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pat)) for pat in alkaloid_patterns)
    if has_alkaloid_pattern:
        return True, "Contains common alkaloid substructure(s)"
    
    return True, "Contains basic nitrogen and heterocyclic rings, within typical molecular weight range"