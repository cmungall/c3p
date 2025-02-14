"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: CHEBI:18366 sphingomyelin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    A sphingomyelin is a phospholipid with a sphingoid base linked to a fatty acid via an amide bond,
    and the terminal hydroxy group of the sphingoid base is esterified to phosphorylcholine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sphingoid base pattern (long aliphatic chain with amino and hydroxy groups)
    sphingoid_pattern = Chem.MolFromSmarts("[NH1X3,NH2X2,NH0X1&r4][CX4H2][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][OX2H]")
    sphingoid_matches = mol.GetSubstructMatches(sphingoid_pattern)
    if len(sphingoid_matches) != 1:
        return False, "No sphingoid base found"

    # Look for fatty acid chain (long carbon chain with carbonyl group)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) != 1:
        return False, "No fatty acid chain found"

    # Check for amide linkage between sphingoid base and fatty acid
    amide_pattern = Chem.MolFromSmarts("[NH1X3,NH2X2,NH0X1&r4]-[CX3](=[OX1])-")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"

    # Look for phosphorylcholine group attached to sphingoid base
    phosphocholine_pattern = Chem.MolFromSmarts("[OX2][CX4][NX4+]([C])(C)[C]")
    phosphocholine_matches = mol.GetSubstructMatches(phosphocholine_pattern)
    if len(phosphocholine_matches) != 1:
        return False, "No phosphorylcholine group found"

    # Check for ester linkage between sphingoid base and phosphorylcholine
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])-")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be sphingomyelin"

    return True, "Contains sphingoid base linked to fatty acid via amide bond, with phosphorylcholine esterified to terminal hydroxy group"