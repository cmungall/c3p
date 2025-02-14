"""
Classifies: CHEBI:36498 galactosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide consists of a galactose head group linked to a ceramide.
    Ceramides have a sphingosine or sphinganine backbone with an amide-linked fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a galactosylceramide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for alpha and beta galactopyranose
    alpha_galactose_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([C@H]([C@H]([C@H](O1)CO)O)O)O)")
    beta_galactose_pattern = Chem.MolFromSmarts("[C@@H]1([C@@H]([C@H]([C@H]([C@H](O1)CO)O)O)O)")
    
    #Check for a galactopyranose moiety.
    if not (mol.HasSubstructMatch(alpha_galactose_pattern) or mol.HasSubstructMatch(beta_galactose_pattern)):
            return False, "No alpha or beta galactose moiety found"

    # SMARTS for sphingosine/sphinganine backbone (including double bond variation)
    ceramide_backbone_pattern_unsaturated = Chem.MolFromSmarts("[C,C@H]([C@H]([C,C@H])O)NC(=O)")
    ceramide_backbone_pattern_saturated = Chem.MolFromSmarts("[C,C@H]([C@H]([C,C@H]([C])([C]))O)NC(=O)")
    
    if not (mol.HasSubstructMatch(ceramide_backbone_pattern_unsaturated) or mol.HasSubstructMatch(ceramide_backbone_pattern_saturated)):
        return False, "No ceramide backbone found"
        
    # Check for glycosidic bond (C-O-C) connecting galactose to ceramide
    glycosidic_bond_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    
    if not glycosidic_matches:
        return False, "No glycosidic bond present"
    
    # Fatty acid chain length check
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "Missing fatty acid chain"


    # Molecular weight check.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
         return False, "Molecular weight too low for galactosylceramide"
    
    return True, "Contains a galactose moiety linked to a ceramide backbone."