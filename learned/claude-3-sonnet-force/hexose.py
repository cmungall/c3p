"""
Classifies: CHEBI:18133 hexose
"""
"""
Classifies: CHEBI:18020 hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is any six-carbon monosaccharide which in its linear form contains 
    either an aldehyde group at position 1 (aldohexose) or a ketone group at position 2 (ketohexose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for 6 carbon atoms and at least 4 oxygen atoms (multiple hydroxyl groups)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count != 6 or o_count < 4:
        return False, "Not a hexose (wrong number of C or O atoms)"
    
    # Check for aldehyde or ketone group
    aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
    ketone_pattern = Chem.MolFromSmarts("C(=O)")
    if not mol.HasSubstructMatch(aldehyde_pattern) and not mol.HasSubstructMatch(ketone_pattern):
        return False, "No aldehyde or ketone group found"
    
    # Check for linear chain or ring structure
    chain_pattern = Chem.MolFromSmarts("[CH2]-[CH]-[CH]-[CH]-[CH]-[CH]")
    pyranose_pattern = Chem.MolFromSmarts("C1OCCCC1")
    furanose_pattern = Chem.MolFromSmarts("C1OCCC1")
    
    if mol.HasSubstructMatch(chain_pattern):
        # Linear hexose, check aldehyde/ketone position
        if mol.HasSubstructMatch(aldehyde_pattern):
            aldehyde_idx = list(mol.GetSubstructMatches(aldehyde_pattern))[0][0]
            if aldehyde_idx != 0:
                return False, "Aldehyde not at position 1 (linear form)"
        else:
            ketone_idx = list(mol.GetSubstructMatches(ketone_pattern))[0][0]
            if ketone_idx != 1:
                return False, "Ketone not at position 2 (linear form)"
    elif mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern):
        # Cyclic hexose, check ring size and arrangement of hydroxyl groups
        pass  # TODO: Implement ring checks
    else:
        return False, "Neither linear nor cyclic hexose structure found"
    
    return True, "Molecule matches the definition of a hexose"