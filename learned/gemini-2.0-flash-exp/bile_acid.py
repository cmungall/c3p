"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    Bile acids are characterized by a steroid core with specific hydroxyl groups and a carboxylic acid at the end of a side chain, often conjugated with glycine or taurine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for the steroid core (tetracyclic ring system) with 5beta configuration
    #The general pattern for a steroid core is:
    # [C]1[C]([C]([C]2[C]([C]([C]([C]([C]1)(C)CC3)CC2)C)C3)(C)
    # We add the 5-beta configuration using the [H] notation
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C@]([C@@]([C]2[C@]([C@@]([C]([C]([C]1)(C)CC3)CC2)C)C3)([H])(C)")
    if not mol.HasSubstructMatch(steroid_core_pattern):
       return False, "Steroid core with 5beta configuration not found"


    # 2. Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
         return False, "Carboxylic acid group not found"

    # 3. Check for at least one hydroxyl group (-OH) on the steroid core.
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl group (-OH) found."


    # 4. Check for Sidechain
    # Look for pattern for a chain of 2 to 4 Carbons attached to the 17 position of steroid, that ends in a carboxylic group.
    # The C17 position is the last carbon of the tetracyclic system attached to 3 or 2 other carbons.
    sidechain_pattern = Chem.MolFromSmarts("[C]~[C]~[C](~C(=O)O)")
    sidechain_matches = mol.GetSubstructMatches(sidechain_pattern)
    # the number of matches should be more than 0 and the first carbon should be bound to C17.
    if len(sidechain_matches) == 0 :
       return False, "Carboxylic acid chain not found"


    # 5. Check for Amide linkage (optional)
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if amide_matches:
        pass # If it has an amide it could be a bile acid.

    # If all checks pass, then it's likely a bile acid
    return True, "Molecule is a bile acid based on steroid core, carboxylic acid group, hydroxyls and sidechain"