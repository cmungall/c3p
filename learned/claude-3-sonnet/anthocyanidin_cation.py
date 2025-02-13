"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: CHEBI:51744 anthocyanidin cation
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    An anthocyanidin cation is an oxygenated derivative of flavylium (2-phenylchromenylium).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for flavylium (2-phenylchromenylium) core
    flavylium_pattern = Chem.MolFromSmarts("[o+]1c2ccccc2cc1")
    if not mol.HasSubstructMatch(flavylium_pattern):
        return False, "Molecule does not contain flavylium core"

    # Check for oxygenated groups (hydroxy, methoxy, etc.)
    has_oxygenated = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen
            has_oxygenated = True
            break
    if not has_oxygenated:
        return False, "Molecule is not oxygenated"

    # Check for positive charge
    if sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() == 1) != 1:
        return False, "Molecule does not have a single positive charge"

    # Check for presence of rings (anthocyanidins are polycyclic)
    if Chem.GetSSSR(mol) == []:
        return False, "Molecule does not contain rings"

    return True, "Molecule contains flavylium core and is oxygenated with a single positive charge"