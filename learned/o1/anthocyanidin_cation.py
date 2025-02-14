"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: anthocyanidin cation
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    An anthocyanidin cation is an aglycone of anthocyanin cation; they are oxygenated derivatives 
    of flavylium (2-phenylchromenylium).
    
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
    
    # Check for positive charge
    charge = rdMolDescriptors.CalcFormalCharge(mol)
    if charge <= 0:
        return False, "Molecule is not a cation"

    # Define SMARTS pattern for flavylium core (2-phenylchromenylium cation)
    # Simplified pattern for flavylium core with positive charge on oxygen
    flavylium_pattern = Chem.MolFromSmarts('c1ccccc1-c2cc3c([o+]cc3)cc2')
    if not mol.HasSubstructMatch(flavylium_pattern):
        return False, "Flavylium core not found"

    # Check for oxygenated substituents (hydroxyl or methoxy groups) on aromatic rings
    oxygenated = False
    oxy_substituent_pattern = Chem.MolFromSmarts('[cH]-[O;H1]')  # Phenolic OH
    methoxy_pattern = Chem.MolFromSmarts('[cH]-O-C')  # Methoxy group

    if mol.HasSubstructMatch(oxy_substituent_pattern):
        oxygenated = True
    elif mol.HasSubstructMatch(methoxy_pattern):
        oxygenated = True

    if not oxygenated:
        return False, "No oxygenated substituents found on aromatic rings"

    # Check for absence of sugar moieties (aglycone)
    # Look for glycosidic bonds: C-O-C between carbons and oxygen
    sugar_pattern = Chem.MolFromSmarts('[C;!R]-O-[C;!R]')
    if mol.HasSubstructMatch(sugar_pattern):
        return False, "Sugar moieties detected (glycosidic bonds present)"

    return True, "Molecule is an anthocyanidin cation"