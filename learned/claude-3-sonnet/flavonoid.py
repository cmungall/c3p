"""
Classifies: CHEBI:47916 flavonoid
"""
Here's a Python program that attempts to classify chemical entities as flavonoids based on their SMILES strings:

"""
Classifies chemical entities of the class CHEBI:26207 flavonoid.

A flavonoid is 'Any member of the 'superclass' flavonoids whose skeleton is based on 1-benzopyran with an aryl substituent at position 2. The term was originally restricted to natural products, but is now also used to describe semi-synthetic and fully synthetic compounds.'
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

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

    # Look for 1-benzopyran core
    benzopyran_pattern = Chem.MolFromSmarts("c1c(oc2ccccc2)cc3ccccc3c1")
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No 1-benzopyran core found"

    # Look for aryl substituent at position 2
    aryl_sub_pattern = Chem.MolFromSmarts("c1c(oc2ccccc2)c(c3ccccc3)cc4ccccc14")
    if not mol.HasSubstructMatch(aryl_sub_pattern):
        return False, "No aryl substituent at position 2 found"

    # Check for common functional groups
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    methoxy_pattern = Chem.MolFromSmarts("[OX2C]")
    glycoside_pattern = Chem.MolFromSmarts("[OX2]C[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)O)O)O)O")

    has_hydroxy = mol.HasSubstructMatch(hydroxy_pattern)
    has_methoxy = mol.HasSubstructMatch(methoxy_pattern)
    has_glycoside = mol.HasSubstructMatch(glycoside_pattern)

    if not (has_hydroxy or has_methoxy or has_glycoside):
        return False, "No common flavonoid functional groups found"

    return True, "Contains 1-benzopyran core with aryl substituent at position 2 and common functional groups"