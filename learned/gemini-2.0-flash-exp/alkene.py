"""
Classifies: CHEBI:32878 alkene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene is an acyclic hydrocarbon with one carbon-carbon double bond
    and general formula CnH2n.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for cycles
    if rdMolDescriptors.CalcNumRings(mol) != 0:
         return False, "Molecule is not acyclic"
    
    # Count double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) != 1:
        return False, f"Molecule has {len(double_bond_matches)} double bond(s). An alkene must have exactly 1."
    
    # Check for only carbon and hydrogen
    c_count = 0
    h_count = 0
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 6:
            c_count += 1
        elif atomic_num == 1:
            h_count += 1
        else:
            return False, "Molecule contains atoms other than C and H"

    # Add explicit hydrogens to molecule
    mol = Chem.AddHs(mol)

    # Count hydrogens again, now with explicit hydrogens
    h_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            h_count += 1

    # Verify formula CnH2n
    if h_count != 2 * c_count:
       return False, f"Formula does not match CnH2n, got C{c_count}H{h_count}"
    
    return True, "Molecule is an acyclic hydrocarbon with one C=C double bond and formula CnH2n"