"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
"""
Classifies: sphingomyelin d18:1
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines if a molecule is a sphingomyelin d18:1 based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a sphingomyelin d18:1, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphocholine head group
    phosphocholine = Chem.MolFromSmarts("[O-]P(=O)(OCC[N+](C)(C)C)O")
    if not mol.HasSubstructMatch(phosphocholine):
        return False, "Missing phosphocholine head group"

    # Check for sphingosine backbone with amide bond
    # [C@H] ensures correct stereochemistry
    # \C=C\ ensures trans double bond
    sphingosine_pattern = Chem.MolFromSmarts("[C@H](O)(/C=C/C)[C@H](COP)NC(=O)")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "Missing characteristic sphingosine backbone structure"

    # Count specific features
    double_bonds = rdMolDescriptors.CalcNumAliphaticDoubleBonds(mol)
    if double_bonds < 1:
        return False, "Missing required double bond"

    # Check for amide group
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "Missing amide bond"

    # Check for two hydroxyl groups (one on sphingosine, one from phosphate)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 1:
        return False, "Missing hydroxyl group on sphingosine"

    # Count carbons in longest chain from NH to end of sphingosine
    # This should be 18 for d18:1
    backbone_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCC/C=C/[C@@H](O)[C@H]")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "Sphingosine backbone must be 18 carbons long (d18:1)"

    # Verify presence of fatty acid chain
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)CCCCC")  # At least 6 carbons
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "Missing fatty acid chain"

    # Check total molecular weight (should be >600 Da for sphingomyelins)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for sphingomyelin"

    # Count key atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if n_count != 2:  # One from amide, one from choline
        return False, "Must have exactly 2 nitrogen atoms"
    if p_count != 1:
        return False, "Must have exactly 1 phosphorus atom"

    return True, "Contains sphingosine d18:1 backbone with phosphocholine head group and fatty acid chain"