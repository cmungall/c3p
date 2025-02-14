"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: CHEBI:35456 flavonoid

A flavonoid is any member of the 'superclass' flavonoids whose skeleton is based on 1-benzopyran
with an aryl substituent at position 2. The term was originally restricted to natural products,
but is now also used to describe semi-synthetic and fully synthetic compounds.
"""

from rdkit import Chem
from rdkit.Chem import rdFMCS

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define flavonoid scaffold pattern
    scaffold_pattern = Chem.MolFromSmarts("[c&r1]1[c&r2]c2c(cc1)oc(c2)-c1ccccc1")

    # Find maximum common substructure with flavonoid scaffold
    mcs = rdFMCS.FindMCS([mol, scaffold_pattern], ringMatchesRingOnly=True)

    # Check if scaffold is present
    if mcs.numAtoms < scaffold_pattern.GetNumAtoms():
        return False, "Flavonoid scaffold not present"

    # Check for aryl substituent at position 2
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and atom.GetTotalNumHs() == 1:
            neighbors = [mol.GetAtomWithIdx(n).GetSymbol() for n in atom.GetNeighbors()]
            if "O" in neighbors and "C" in neighbors:
                aryl_sub = True
                break
    else:
        return False, "No aryl substituent at position 2"

    return True, "Contains flavonoid scaffold with aryl substituent at position 2"