"""
Classifies: CHEBI:18310 alkane
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkane(smiles: str):
    """
    Determines if a molecule is an alkane based on its SMILES string.
    An alkane is an acyclic, saturated hydrocarbon with the general formula CnH2n+2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkane, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check if only carbon and hydrogen atoms are present
    for atom in mol.GetAtoms():
      atomic_num = atom.GetAtomicNum()
      if atomic_num != 1 and atomic_num != 6:
        return False, "Molecule contains non-carbon and non-hydrogen atoms"

    # 2. Check if all bonds are single bonds
    for bond in mol.GetBonds():
      if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
        return False, "Molecule contains double or triple bonds"
    
    # 3. Check if the molecule is acyclic (no rings)
    if rdMolDescriptors.CalcNumRings(mol) > 0:
        return False, "Molecule contains rings"

    #4 No need to check the formula, as we are implicitly checking for saturation
    # and presence of only C and H

    return True, "Molecule is an acyclic, saturated hydrocarbon (alkane)"