"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for glycerol backbone with proper chiral center
    glycerol_pattern = Chem.MolFromSmarts("O[C@@H](COC(=O)[C,H])COP(=O)(O)OCCN")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone with chiral center and phosphate-ethanolamine linkage not found"

    # Exclude common motifs of similar but non-matching phospholipids
    wrong_substructs = [
        Chem.MolFromSmarts("O=C(OCC[N+](C)(C)C)P(=O)(O)OCC"),
        Chem.MolFromSmarts("OC[C@H](NP(=O)(OC)O)CO")
    ]
    for sub in wrong_substructs:
        if mol.HasSubstructMatch(sub):
            return False, "Found structural elements common to other phospholipids like PC or PS"
    
    # Ensure 1-O-acyl group is linked in correct position
    acyl_pattern = Chem.MolFromSmarts("OC(=O)[C,H]")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No 1-O-acyl ester linkage found"
    
    return True, "Molecule fits all structural criteria for 1-acyl-sn-glycero-3-phosphoethanolamine"