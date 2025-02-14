"""
Classifies: CHEBI:60834 beta-carbolines
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.
    A beta-carboline is a pyrido[3,4-b]indole or its hydrogenated derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-carboline, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core beta-carboline structure using SMARTS. Allow for single or double bonds.
    # The N in the indole should be connected to a ring carbon to be a true indole
    core_pattern1 = Chem.MolFromSmarts('c1cc2[nH]c3C[C,c]c(C[C,c]2c1)[C,c]3') #allow for fully aromatic
    core_pattern2 = Chem.MolFromSmarts('c1cc2[nH]c3CCN[C,c](C[C,c]2c1)[C,c]3') #allow for double bond reduction
    core_pattern3 = Chem.MolFromSmarts('c1cc2[nH]c3C[C,c]N[C,c](C[C,c]2c1)[C,c]3') #allow for both reductions
    core_pattern4 = Chem.MolFromSmarts('c1cc2nc3C[C,c]c(C[C,c]2c1)[C,c]3') #allow for fused pyridine
    core_pattern5 = Chem.MolFromSmarts('c1cc2nc3CCN[C,c](C[C,c]2c1)[C,c]3')# allow for fused pyridine and reduction
    core_pattern6 = Chem.MolFromSmarts('c1cc2nc3C[C,c]N[C,c](C[C,c]2c1)[C,c]3')# allow for fused pyridine and both reductions
    core_pattern7 = Chem.MolFromSmarts('c1cc2[nH]c3[CH2][CH2][C,c](C[C,c]2c1)[C,c]3') #allow for dihydropyrido
    core_pattern8 = Chem.MolFromSmarts('c1cc2nc3[CH2][CH2][C,c](C[C,c]2c1)[C,c]3') #allow for dihydropyrido, fused pyridine
    core_pattern9 = Chem.MolFromSmarts('C1CC2=C3C(=CC=C1)NC1=CC=CC=C3N21') #another common form


    if (mol.HasSubstructMatch(core_pattern1) or
        mol.HasSubstructMatch(core_pattern2) or
        mol.HasSubstructMatch(core_pattern3) or
        mol.HasSubstructMatch(core_pattern4) or
        mol.HasSubstructMatch(core_pattern5) or
        mol.HasSubstructMatch(core_pattern6) or
        mol.HasSubstructMatch(core_pattern7) or
        mol.HasSubstructMatch(core_pattern8) or
        mol.HasSubstructMatch(core_pattern9)):
      return True, "Contains the beta-carboline core structure"

    return False, "Does not contain the beta-carboline core structure"