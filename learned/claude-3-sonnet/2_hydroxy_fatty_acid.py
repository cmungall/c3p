"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: CHEBI:51865 2-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid is any fatty acid with a hydroxy functional group in the 2-position.

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

    # Look for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for hydroxy group in the 2-position relative to the acid group
    hydroxy_pattern = Chem.MolFromSmarts("[C;H2]([C;H1]([OH]))[C;H2]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if not hydroxy_matches:
        return False, "No 2-hydroxy group found"

    # Check for long carbon chain (fatty acid)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "No fatty acid chain found"

    # Check for stereochemistry
    atoms = [mol.GetAtomWithIdx(idx) for idx in hydroxy_matches[0]]
    if any(atom.HasProp("_CIPCode") for atom in atoms):
        stereo_pattern = Chem.MolFromSmarts("[C@H](O)[C@H](C)C")
        if not mol.HasSubstructMatch(stereo_pattern):
            return False, "Incorrect stereochemistry for 2-hydroxy group"

    # Additional checks
    mol_wt = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, "Molecular weight too low for a fatty acid"

    return True, "Contains a 2-hydroxy group in a fatty acid chain"