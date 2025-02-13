"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: CHEBI:16540 sphingomyelin

A sphingomyelin is a phospholipid where the amino group of a sphingoid base is in amide linkage with
one of several fatty acids, while the terminal hydroxy group of the sphingoid base is esterified to
phosphorylcholine.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.

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
    
    # Look for sphingoid base pattern (long hydrocarbon chain with amino and hydroxy groups)
    sphingoid_pattern = Chem.MolFromSmarts("C[C@@H](O)[C@@H](N)CCCCCCCCCCCCC")
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid base found"
    
    # Look for fatty acid chain (long hydrocarbon chain attached to amino group via amide)
    fatty_acid_pattern = Chem.MolFromSmarts("CCC(=O)N[C@@H](CCCCCCCCCCCCC)[C@@H](O)CCCCCCCCCCCCC")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No fatty acid chain attached via amide linkage"
    
    # Look for phosphorylcholine group (ester-linked to terminal hydroxy of sphingoid base)
    phosphocholine_pattern = Chem.MolFromSmarts("COP(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphorylcholine group found"
    
    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be sphingomyelin"
    
    # Check molecular weight - sphingomyelins typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for sphingomyelin"
    
    # Count carbons, nitrogens, oxygens and phosphorus
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if c_count < 30:
        return False, "Too few carbons for sphingomyelin"
    if n_count != 2:
        return False, "Must have exactly 2 nitrogens"
    if o_count < 5:
        return False, "Too few oxygens for sphingomyelin"
    if p_count != 1:
        return False, "Must have exactly 1 phosphorus"
    
    return True, "Contains a sphingoid base with fatty acid chain attached via amide, and phosphorylcholine ester"