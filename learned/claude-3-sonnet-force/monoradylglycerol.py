"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: CHEBI:36225 monoradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol is a glycerol bearing a single acyl, alkyl or alk-1-enyl substituent at an unspecified position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for ester (-O-C(=O)-), ether (-O-), or alkene (-C=C-) group
    substituent_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])|[OX2][CX4]|[CX3]=[CX3]")
    substituent_matches = mol.GetSubstructMatches(substituent_pattern)
    if len(substituent_matches) != 1:
        return False, f"Found {len(substituent_matches)} substituents, need exactly 1"
    
    # Count rotatable bonds to verify chain length
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Substituent chain too short for monoradylglycerol"
    
    # Check molecular weight - monoradylglycerols typically >250 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low for monoradylglycerol"
    
    # Count carbons, oxygens, and double bonds
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_double_bonds = rdMolDescriptors.CalcNumAromaticRings(mol) + rdMolDescriptors.CalcNumRotatableBonds(mol) - n_rotatable
    
    if c_count < 10:
        return False, "Too few carbons for monoradylglycerol"
    if o_count != 4:
        return False, "Must have exactly 4 oxygens (3 in glycerol, 1 in substituent)"
    if n_double_bonds > 1:
        return False, "Too many double bonds for monoradylglycerol"
        
    return True, "Contains glycerol backbone with a single acyl, alkyl or alk-1-enyl substituent"