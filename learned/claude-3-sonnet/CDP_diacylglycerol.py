"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
"""
Classifies: CDP-diacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a CDP-diacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for cytosine base connected to ribose sugar
    cytidine_patterns = [
        # Various representations of cytidine (cytosine + ribose)
        Chem.MolFromSmarts("[CH2]1[CH]([CH]([OH])[CH]([OH])O1)N2C=CC(=NC2=O)N"),
        Chem.MolFromSmarts("[CH2]1[CH]([CH]([OH])[CH]([OH])O1)n2ccc(N)nc2=O"),
        Chem.MolFromSmarts("OC[C@H]1O[C@@H](N2C=CC(=NC2=O)N)C(O)[C@H]1O"),
        Chem.MolFromSmarts("OC[C@H]1O[C@H](n2ccc(N)nc2=O)[C@H](O)[C@@H]1O")
    ]
    
    has_cytidine = any(mol.HasSubstructMatch(pattern) for pattern in cytidine_patterns if pattern is not None)
    if not has_cytidine:
        return False, "No cytidine moiety found"

    # Check for diphosphate bridge
    diphosphate_patterns = [
        Chem.MolFromSmarts("OP(=O)(O)OP(=O)(O)O"),
        Chem.MolFromSmarts("[OX2]P(=O)([OX2])OP(=O)([OX2])[OX2]"),
        Chem.MolFromSmarts("P(O)(=O)OP(O)(=O)O")
    ]
    has_diphosphate = any(mol.HasSubstructMatch(pattern) for pattern in diphosphate_patterns if pattern is not None)
    if not has_diphosphate:
        return False, "No diphosphate bridge found"

    # Check for glycerol backbone with two ester groups
    glycerol_patterns = [
        Chem.MolFromSmarts("[CH2X4]([OX2]C(=O))[CHX4]([OX2]C(=O))[CH2X4]O"),
        Chem.MolFromSmarts("[CH2]([OX2]C(=O))[CH]([OX2]C(=O))[CH2]OP"),
        Chem.MolFromSmarts("OCC(COC(=O))OC(=O)"),
        Chem.MolFromSmarts("OC[C@H](OC(=O))COC(=O)")
    ]
    has_glycerol = any(mol.HasSubstructMatch(pattern) for pattern in glycerol_patterns if pattern is not None)
    if not has_glycerol:
        return False, "No glycerol backbone with two ester groups found"

    # Count ester groups (should be exactly 2)
    ester_pattern = Chem.MolFromSmarts("[OX2]C(=O)[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check for long carbon chains (fatty acids)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, "Missing long carbon chains for fatty acids"

    # Verify complete connectivity
    # The molecule should have all components connected:
    # cytidine-diphosphate-glycerol-fatty_acids
    if len(Chem.GetMolFrags(mol)) > 1:
        return False, "Disconnected fragments found"

    # Count key atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)

    # Basic element count checks
    if c_count < 20:
        return False, "Too few carbons for CDP-diacylglycerol"
    if o_count < 12:
        return False, "Too few oxygens for CDP-diacylglycerol"
    if n_count != 3:
        return False, "Should have exactly 3 nitrogens (cytosine)"
    if p_count != 2:
        return False, "Should have exactly 2 phosphorus atoms"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800:  # CDP-diacylglycerols are typically >800 Da
        return False, "Molecular weight too low for CDP-diacylglycerol"

    return True, "Contains cytidine diphosphate group connected to glycerol backbone with two fatty acid chains"