"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: Phenylpropanoid
"""

from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    Phenylpropanoids are organic aromatic compounds with a structure based on a phenylpropane skeleton (C6-C3).
    This includes subclasses such as flavonoids, coumarins, lignans, and anthocyanins.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Phenylpropanoid general pattern: benzene ring connected to a three-carbon chain (C6-C3 unit)
    # Allow for variations, including unsaturation and substitution in the propanoid chain
    phenylpropanoid_pattern = Chem.MolFromSmarts("c1ccccc1CCC")
    phenylpropenoid_pattern = Chem.MolFromSmarts("c1ccccc1CC=C")
    cinnamic_acid_pattern = Chem.MolFromSmarts("c1ccccc1C=CC(=O)O")
    cinnamyl_alcohol_pattern = Chem.MolFromSmarts("c1ccccc1C=CCO")

    # Flavonoid core structure (C6-C3-C6 skeleton with heterocyclic ring)
    flavonoid_pattern = Chem.MolFromSmarts("c1cc(-c2ccc3occ(=O)c(-c4ccccc4)c3c2)ccc1")

    # Coumarin core structure
    coumarin_pattern = Chem.MolFromSmarts("O=C1C=CC2=CC=CC=C2O1")

    # Lignan structure (two phenylpropane units linked)
    lignan_pattern = Chem.MolFromSmarts("c1ccccc1C[C@@H](C)c2ccc(O)cc2")  # Simplified pattern

    # Anthocyanin core structure (flavylium ion)
    anthocyanin_pattern = Chem.MolFromSmarts("[O+]c1cc2ccccc2oc1")

    # Stilbenoid core structure
    stilbenoid_pattern = Chem.MolFromSmarts("c1ccc(cc1)/C=C/c2ccccc2")

    # List of patterns to check with corresponding explanations
    patterns = [
        ("phenylpropane skeleton", phenylpropanoid_pattern),
        ("phenylpropenoid skeleton", phenylpropenoid_pattern),
        ("cinnamic acid structure", cinnamic_acid_pattern),
        ("cinnamyl alcohol structure", cinnamyl_alcohol_pattern),
        ("flavonoid core structure", flavonoid_pattern),
        ("coumarin core structure", coumarin_pattern),
        ("lignan structure", lignan_pattern),
        ("anthocyanin core structure", anthocyanin_pattern),
        ("stilbenoid core structure", stilbenoid_pattern),
    ]

    for (description, pattern) in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains {description}"

    return False, "Does not contain phenylpropanoid skeleton"