"""
Classifies: CHEBI:48039 dihydroflavonols
"""
"""
Classifies: CHEBI:15864 dihydroflavonol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol is a hydroxyflavanone with a hydroxy group at position 3
    of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic flavanone core structure with 3-OH
    # [#6]1-[#6]-[#6](=[#8])-c2c(O1)cccc2
    # Adding 3-OH position requirement
    dihydroflavonol_pattern = Chem.MolFromSmarts(
        '[#6]1-[#6]([OH1])-[#6](=[#8])-c2c(O1)cccc2'
    )
    
    if not mol.HasSubstructMatch(dihydroflavonol_pattern):
        return False, "Missing dihydroflavonol core structure with 3-OH group"

    # Must have phenyl ring attached at position 2
    phenyl_pattern = Chem.MolFromSmarts('O1[CH1][c]2[c,n]cccc2')
    if not mol.HasSubstructMatch(phenyl_pattern):
        return False, "Missing 2-phenyl substituent"

    # Count carbons and oxygens to ensure basic composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 15:
        return False, "Too few carbons for dihydroflavonol structure"
    if o_count < 4:
        return False, "Too few oxygens for dihydroflavonol structure"

    # Check for ketone group in the correct position
    ketone_pattern = Chem.MolFromSmarts('[#6]-[#6](=[#8])-[#6]')
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "Missing ketone group"

    # Verify the presence of the chromanone system
    chromanone_pattern = Chem.MolFromSmarts('O1[#6][#6][#6](=[#8])c2c1cccc2')
    if not mol.HasSubstructMatch(chromanone_pattern):
        return False, "Missing chromanone system"

    return True, "Contains dihydroflavonol core structure with 3-OH group"