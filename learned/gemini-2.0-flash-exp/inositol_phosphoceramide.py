"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    An inositol phosphoceramide is a phosphosphingolipid with an inositol and ceramide
    linked via a phosphodiester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an inositol phosphoceramide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Ceramide Identification
    # Look for amide bond (C(=O)N), 2 hydroxyls and a carbon chain (sphingosine)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond found, not a ceramide"

    # Identify the sphingosine backbone
    sphingosine_pattern = Chem.MolFromSmarts("[CX4](O)[CX4](O)")
    sphingosine_matches = mol.GetSubstructMatches(sphingosine_pattern)
    if not sphingosine_matches:
      return False, "No sphingosine backbone found"
    
    # 2. Phosphate group
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"

    # 3. Inositol Identification
    # Look for the inositol ring with 6 hydroxyls
    inositol_pattern = Chem.MolFromSmarts("[C]1([O])[C]([O])[C]([O])[C]([O])[C]([O])[C]1[O]")
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_matches:
        return False, "No inositol ring found"
    
    # 4. Verify linkage
    # Check if phosphate links to inositol
    # And check if phosphate also links to the ceramide (via C-O-P)
    inositol_phos_pattern = Chem.MolFromSmarts("[C]1([O])[C]([O])[C]([O])[C]([O])[C]([O])[C]1[O][OP](=O)(O)[OX2][CX4]")
    inositol_phos_matches = mol.GetSubstructMatches(inositol_phos_pattern)
    if not inositol_phos_matches:
        return False, "Phosphate is not directly linked to inositol"

    ceramide_phos_pattern = Chem.MolFromSmarts("[NX3][CX4]([OX2])[CX4][OX2][P](=[OX1])([OX2])")
    ceramide_phos_matches = mol.GetSubstructMatches(ceramide_phos_pattern)
    if not ceramide_phos_matches:
      return False, "Phosphate is not directly linked to ceramide"

    #5. Optional: check for mannose
    mannose_pattern = Chem.MolFromSmarts("OC[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)O1")
    mannose_matches = mol.GetSubstructMatches(mannose_pattern)
    
    # If all checks pass, return True
    return True, "Structure contains ceramide, inositol, and phosphate group with correct linkage"