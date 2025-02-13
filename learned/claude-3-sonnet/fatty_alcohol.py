"""
Classifies: CHEBI:24026 fatty alcohol
"""
"""
Classifies: CHEBI:36195 fatty alcohol

A fatty alcohol is an aliphatic alcohol consisting of a chain of 3 to greater than
27 carbon atoms. Fatty alcohols may be saturated or unsaturated and may be branched
or unbranched.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count alcohol groups (-OH)
    alcohol_pattern = Chem.MolFromSmarts("[OX1H]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    if len(alcohol_matches) != 1:
        return False, f"Found {len(alcohol_matches)} alcohol groups, need exactly 1"

    # Count carbon chain length
    carbon_chain = max(len(AllChem.GenerateDepictionMatching3DStructure(mol, Chem.MolFromSmarts("[C]~[C]"), useHs=False)) + 1
    if carbon_chain < 3 or carbon_chain > 27:
        return False, f"Carbon chain length ({carbon_chain}) outside of allowed range (3-27)"

    # Check for aliphatic (no cyclic structures)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains cyclic structures, must be aliphatic"

    # Count unsaturations (double bonds)
    num_double_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol) - rdMolDescriptors.CalcNumRotatableBonds(Chem.RemoveHs(mol))
    if num_double_bonds > carbon_chain - 3:
        return False, "Too many unsaturations for carbon chain length"

    return True, "Aliphatic alcohol with 3-27 carbon chain length"