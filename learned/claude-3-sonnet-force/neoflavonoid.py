"""
Classifies: CHEBI:71971 neoflavonoid
"""
"""
Classifies: CHEBI:35706 neoflavonoid

A neoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 4. 
The term was originally restricted to natural products, but is now also used to describe 
semi-synthetic and fully synthetic compounds.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for benzopyran core
    benzopyran_pattern = Chem.MolFromSmarts("c1c(oc2ccccc2)cccc1")
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No benzopyran core found"

    # Check for aryl substituent at position 4
    aryl_pattern = Chem.MolFromSmarts("[cX3](-[cX3])-[cX3]")
    aryl_matches = mol.GetSubstructMatches(aryl_pattern)
    if not aryl_matches:
        return False, "No aryl substituent at position 4"

    # Check if aryl substituent is attached to benzopyran at position 4
    for match in aryl_matches:
        atom_idx = match[0]
        env = Chem.FindAtomEnvironmentOfRadiusN(mol, 4, atom_idx)
        if env.getIsConstrained():
            # Check if benzopyran core is part of environment
            if benzopyran_pattern.HasSubstructMatch(env.Abbreviation):
                return True, "Molecule has a benzopyran core with an aryl substituent at position 4"

    return False, "Aryl substituent not attached to benzopyran at position 4"