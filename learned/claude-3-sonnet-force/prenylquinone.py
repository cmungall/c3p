"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: CHEBI:51333 prenylquinone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone is a quinone substituted by a polyprenyl-derived side-chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of quinone substructure
    quinone_pattern = Chem.MolFromSmarts("[#6]1=,:[#6](~[#6]=,:[#6](~[#6]1~[#8])~[#8])~[#8]")
    if not mol.HasSubstructMatch(quinone_pattern):
        return False, "No quinone substructure found"
    
    # Check for presence of prenyl side chain (at least one)
    prenyl_pattern = Chem.MolFromSmarts("[#6]-[#6]=,:[#6]-[#6](-[#6])(-[#6])-[#6]")
    if not mol.HasSubstructMatch(prenyl_pattern):
        return False, "No prenyl side chain found"
    
    # Check for long carbon chains (at least 6 carbons in a row)
    long_chain_pattern = Chem.MolFromSmarts("[#6]-[#6]-[#6]-[#6]-[#6]-[#6]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long carbon chains found"
    
    # Check molecular weight - prenylquinones typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for prenylquinone"
    
    return True, "Contains quinone substructure and prenyl side chain(s)"