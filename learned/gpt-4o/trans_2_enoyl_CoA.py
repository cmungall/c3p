"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
from rdkit import Chem

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trans-2-enoyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for general adenine structure (part of CoA)
    coa_adenine_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)N")
    if not mol.HasSubstructMatch(coa_adenine_pattern):
        return False, "Not a CoA molecule, adenine not found"

    # Check for the generalized CoA structure (adenine, ribose, phosphate chains)
    # Using a more general smart pattern for CoA core that matches essential functional groups with more variability
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(O)=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Full CoA structure not confirmed"

    # Ensure a thioester linkage is present, allowing for variability in aliphatic chain
    thiol_ester_pattern = Chem.MolFromSmarts("C(=O)SC")
    if not mol.HasSubstructMatch(thiol_ester_pattern):
        return False, "Missing thioester linkage"

    # Check specifically for the trans double bond at the second position
    # Broad trans configuration capturing both E/Z nomenclature
    trans_double_bond_pattern = Chem.MolFromSmarts("C/C=C\\C")
    double_bond_matches = mol.GetSubstructMatches(trans_double_bond_pattern)
    if len(double_bond_matches) == 0:
        return False, "No 2,3-trans double bond found"

    return True, "Molecule fits the trans-2-enoyl-CoA structure criteria"