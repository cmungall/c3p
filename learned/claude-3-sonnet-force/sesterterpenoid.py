"""
Classifies: CHEBI:26660 sesterterpenoid
"""
"""
Classifies: CHEBI:33551 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    A sesterterpenoid is a terpenoid with a C25 backbone, derived from a sesterterpene.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 25:
        return False, f"Molecule has {c_count} carbon atoms, should have 25"
    
    # Check for terpenoid-like structure
    try:
        Chem.AddHs(mol)
        Chem.EmbedMolecule(mol)
        Chem.AssignStereochemistry(mol)
    except ValueError:
        return False, "Could not assign stereochemistry"
    
    # Check for rings and double bonds
    n_rings = mol.GetRingInfo().NumRings()
    n_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    
    # Sesterterpenoids typically have 4-7 rings and 4-8 double bonds
    if n_rings < 4 or n_rings > 7:
        return False, f"Unexpected number of rings ({n_rings}), should be 4-7"
    if n_double_bonds < 4 or n_double_bonds > 8:
        return False, f"Unexpected number of double bonds ({n_double_bonds}), should be 4-8"
    
    # Check molecular weight - sesterterpenoids typically 300-500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 500:
        return False, f"Molecular weight ({mol_wt:.2f} Da) out of typical sesterterpenoid range (300-500 Da)"
    
    return True, "Molecule matches structural characteristics of a sesterterpenoid"