"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: CHEBI:35510 acyl-CoA
An acyl-CoA is a thioester that results from the formal condensation of the thiol group of coenzyme A
with the carboxy group of any carboxylic acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA contains a thioester linkage between a carboxylic acid and the thiol group of Coenzyme A.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for thioester pattern (-C(=O)-S-)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Look for CoA backbone pattern
    coa_pattern = Chem.MolFromSmarts("C1(OP(OP(OCC(C(NC(CCNC(=O)CCNC(=O)[C@H](O)[C@@H](C)[C@@H](COP(O)(O)=O)OP(O)(O)=O)[n+]1cnc2c(N)ncnc2N)=O)=O)C)C)(O)O)[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA backbone found"
    
    # Check that thioester links carboxylic acid to CoA
    for match in thioester_matches:
        acid_idx = mol.GetAtomWithIdx(match[0]).GetNeighbors()[0].GetIdx()
        coa_idx = mol.GetAtomWithIdx(match[1]).GetNeighbors()[0].GetIdx()
        
        acid_mol = Chem.MolFromSmiles(Chem.MolFragmentToSmiles(mol, [acid_idx], bondsToUse=[]))
        coa_mol = Chem.MolFromSmiles(Chem.MolFragmentToSmiles(mol, [coa_idx], bondsToUse=[]))
        
        acid_smarts = Chem.MolToSmarts(acid_mol)
        coa_smarts = Chem.MolToSmarts(coa_mol)
        
        if acid_smarts != "C(=O)O" or not coa_mol.HasSubstructMatch(coa_pattern):
            return False, "Thioester does not link carboxylic acid and CoA"

    return True, "Contains thioester linkage between carboxylic acid and CoA thiol group"