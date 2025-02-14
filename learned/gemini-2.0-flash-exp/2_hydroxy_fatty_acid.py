"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid has a carboxylic acid group and a hydroxyl group
    on the carbon alpha to the carboxyl carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
         return False, "No carboxylic acid group found"

    # Check for 2-hydroxy group (hydroxy on the alpha carbon)
    hydroxy_pattern = Chem.MolFromSmarts("[CHX4]([OH])[C](=O)O")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No hydroxyl group at the alpha position found"

    # Check for long chain - count non-acid carbon atoms > 3
    carbon_count = 0
    for atom in mol.GetAtoms():
      if atom.GetAtomicNum() == 6:
        is_acid_carbon = False
        for bond in atom.GetBonds():
           other_atom = bond.GetOtherAtom(atom)
           if other_atom.GetAtomicNum() == 8 and other_atom.GetIdx() != atom.GetIdx():
             is_acid_carbon = True
             break
        if not is_acid_carbon:
          carbon_count += 1
    if carbon_count < 3:
         return False, "Carbon chain too short for a fatty acid"

    # Additional check to remove edge cases like simple molecules.
    if rdMolDescriptors.CalcNumHeavyAtoms(mol) < 5:
        return False, "Molecule too small to be a fatty acid"

    return True, "Contains a carboxylic acid group with a hydroxyl group on the alpha carbon and a carbon chain"