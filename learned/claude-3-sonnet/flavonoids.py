"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: CHEBI:25107 flavonoids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are characterized by a C15 or C16 skeleton with a phenyl-substituted 1-phenylpropane structure,
    or a condensed C6-C3 lignan precursor.

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
    
    # Define flavonoid core structure patterns
    flavonoid_core_patterns = [
        Chem.MolFromSmarts("[c1c(cc2c(c1)OCc3c(c2)cccc3)O]"),   # Flavone
        Chem.MolFromSmarts("[c1c(cc2c(c1)OCC3=CC(=O)Oc4c3cccc4)O]"),  # Flavonol
        Chem.MolFromSmarts("[c1c(cc2c(c1)OCC3=CC(=O)c4c(O)cccc4O3)O]"),  # Flavanone
        Chem.MolFromSmarts("[c1c(c2c(cc1O)OCC3=CC(=O)c4c(O)cccc4O3)O]"),  # Flavan-3-ol
        Chem.MolFromSmarts("[c1c(c2c(cc1O)OCC3=C(O)c4c(cccc4O3)O)O]"),  # Flavan-4-ol
        Chem.MolFromSmarts("[c1c(c2c(cc1O)OCC3=C(O)c4c(cc(O)cc4O3)O)O]"),  # Flavan-3,4-diol
        Chem.MolFromSmarts("[c1c(c2c(cc1O)OCC3=CC(=O)c4c(ccc(O)c4O3)O)O]"),  # Isoflavone
        Chem.MolFromSmarts("[c1c(c2c(cc1O)Oc3c(c2)cccc3)O]"),  # Chalcone
        Chem.MolFromSmarts("[c1c(c2c(cc1O)Oc3c(c2)c(O)ccc3)O]"),  # Aurone
        Chem.MolFromSmarts("[c1c(c2c(cc1O)Oc3c(c2)cccc3O)O]"),  # Dihydrochalcone
        Chem.MolFromSmarts("[c1c(c2c(cc1O)Oc3c(c2)cccc3)CO]"),  # Neoflavonoid
    ]
    
    # Check for flavonoid core structure
    for pattern in flavonoid_core_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Molecule contains a flavonoid core structure"
    
    # Check for lignan precursor
    lignan_pattern = Chem.MolFromSmarts("[c1c(ccc2c1Oc3ccccc3O2)O]")
    if mol.HasSubstructMatch(lignan_pattern):
        return True, "Molecule contains a condensed C6-C3 lignan precursor"
    
    return False, "No flavonoid or lignan precursor structure found"