"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    An N-acylsphingosine has a sphingosine backbone with a fatty acyl group attached to the nitrogen via an amide bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphingosine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for sphingosine backbone, including stereochemistry
    sphingosine_pattern = Chem.MolFromSmarts("[C@H]([C@@H](O)[H])([H])(NC)[C@@H](O)[C@@H]([H])\C=C\CCCCCCCCCCCCC")
    if not mol.HasSubstructMatch(sphingosine_pattern):
       return False, "No sphingosine backbone found"

    # Define SMARTS for the amide bond
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)

    if len(amide_matches) !=1:
        return False, f"Found {len(amide_matches)} amide bonds, need exactly 1"

    # Check for long carbon chain on the acyl side
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4][CX4]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) == 0:
        return False, f"Missing fatty acid chain"

    # Check molecular weight - N-acylsphingosines typically >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
      return False, "Molecular weight too low for N-acylsphingosine"

    # Check for a long carbon chain in sphingosine (12 carbon in the pattern) and also acyl chain
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 18:
       return False, "Too few carbons for N-acylsphingosine"

    # Verify there are 2 OH groups
    oh_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and len(atom.GetNeighbors()) == 1)
    if oh_count < 2:
        return False, "Must have at least 2 hydroxyl groups"
    
    return True, "Contains a sphingosine backbone with a fatty acyl group attached via an amide bond"