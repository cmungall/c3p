"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: CHEBI:60367 nucleoside 5'-phosphate

A nucleoside 5'-phosphate is a ribosyl or deoxyribosyl derivative of a pyrimidine or purine base
in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for ribose/deoxyribose backbone
    ribose_pattern = Chem.MolFromSmarts("[OX2]C1[CH]([OX2])[CH]([OX2])[CH]([OX2])[CH]1")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose/deoxyribose backbone found"
    
    # Look for pyrimidine or purine base
    base_patterns = [Chem.MolFromSmarts(sma) for sma in ["c1ncnc2n1cncn12", "c1cnc2n1cncn12",
                                                          "c1ncnc2[nH]1ccn12", "c1cnc2[nH]1ccn12"]]
    base_match = any(mol.HasSubstructMatch(pat) for pat in base_patterns)
    if not base_match:
        return False, "No pyrimidine or purine base found"
    
    # Look for phosphate group(s) attached to C-5 of ribose
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    ribose_atoms = mol.GetSubstructMatches(ribose_pattern)[0]
    c5_idx = ribose_atoms[4]
    phosphate_attached = any(c5_idx in match for match in phosphate_matches)
    
    if not phosphate_attached:
        return False, "No phosphate group attached to C-5 of ribose"
    
    return True, "Contains ribose/deoxyribose backbone with pyrimidine/purine base and phosphate group on C-5"