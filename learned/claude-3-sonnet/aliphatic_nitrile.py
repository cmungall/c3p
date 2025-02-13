"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: CHEBI:46907 aliphatic nitrile
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile is any nitrile derived from an aliphatic compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for nitrile group (-C≡N)
    nitrile_pattern = Chem.MolFromSmarts("[C;X2]#[N;X1]")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)
    if not nitrile_matches:
        return False, "No nitrile group found"

    # Check if the molecule is aliphatic (no aromatic rings)
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings > 0:
        return False, "Molecule contains aromatic rings, not aliphatic"

    # Check if the molecule contains only C, H, N, and optionally O, S, P, F, Cl, Br, I
    allowed_atoms = set([6, 1, 7, 8, 16, 15, 9, 17, 35, 53])
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Molecule contains disallowed atom {atom.GetSymbol()}"

    return True, "Molecule is an aliphatic nitrile"