"""
Classifies: CHEBI:28892 ganglioside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    A ganglioside is a glycosphingolipid containing one or more sialic acids 
    linked to an oligosaccharide chain attached to a ceramide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight - gangliosides are large molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800:
        return False, "Molecular weight too low for ganglioside"

    # Look for ceramide backbone
    # Sphingosine pattern (long chain base with OH and NH)
    sphingosine_pattern = Chem.MolFromSmarts("[CH2]-[CH]-[CH](-[OH])-[CH](-[NH]-[C](=O))-[CH2]-O")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No ceramide backbone found"

    # Look for sialic acid (Neu5Ac) pattern
    # Core structure of sialic acid with carboxyl group
    sialic_pattern = Chem.MolFromSmarts("[C](=O)[O;H,-]-[C]1-[C]-[C](-[NH]-[C](=O)-[CH3])-[C]-[C](-[OH])-[C]1")
    sialic_matches = mol.GetSubstructMatches(sialic_pattern)
    if len(sialic_matches) < 1:
        return False, "No sialic acid (Neu5Ac) found"

    # Look for oligosaccharide pattern
    # Basic sugar rings with OH groups
    sugar_pattern = Chem.MolFromSmarts("[C]1-[O]-[C]-[C]-[C]-[C]-1")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) < 3:
        return False, "Insufficient oligosaccharide chain"

    # Check for glycosidic linkages
    glycosidic_pattern = Chem.MolFromSmarts("[C]-[O]-[C]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if len(glycosidic_matches) < 4:
        return False, "Insufficient glycosidic linkages"

    # Count key atoms to verify overall composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 30:
        return False, "Too few carbons for ganglioside structure"
    if o_count < 15:
        return False, "Too few oxygens for ganglioside structure"
    if n_count < 2:
        return False, "Too few nitrogens for ganglioside structure"

    # Verify presence of long chain fatty acid
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "Missing long chain fatty acid"

    return True, "Contains ceramide backbone, sialic acid(s), and oligosaccharide chain with correct linkages"