"""
Classifies: CHEBI:67194 cannabinoid
"""
from rdkit import Chem

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    This function uses specific substructure searches rather than a scoring system
    to identify cannabinoid characteristics.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for dibenzopyran core (classical cannabinoids) - stricter pattern
    dibenzopyran_pattern = Chem.MolFromSmarts("c1cc2c3c(cc1)Oc1ccccc1C3CC2")
    if mol.HasSubstructMatch(dibenzopyran_pattern):
        # Check for at least one phenolic -OH
        phenolic_oh_pattern = Chem.MolFromSmarts("[c]O[H]")
        if mol.HasSubstructMatch(phenolic_oh_pattern):
          return True, "Contains classical cannabinoid dibenzopyran core and a phenolic hydroxyl."

    # 2. Check for modified dibenzopyran core
    modified_dibenzopyran_pattern = Chem.MolFromSmarts("c1cc2c(cc1[C,c])Oc1ccccc1C2")
    if mol.HasSubstructMatch(modified_dibenzopyran_pattern):
       # Check for at least one phenolic -OH
      phenolic_oh_pattern = Chem.MolFromSmarts("[c]O[H]")
      if mol.HasSubstructMatch(phenolic_oh_pattern):
        return True, "Contains modified dibenzopyran core and a phenolic hydroxyl."

    # 3. Check for a cyclic terpene core with phenolic OH. For CBD like compounds.
    cyclic_terpene_pattern = Chem.MolFromSmarts("[C]1(C)C=C([C,c])[C,c]CC1") # Simplified, look for monoterpene core
    if mol.HasSubstructMatch(cyclic_terpene_pattern):
        # Check for at least one phenolic -OH
      phenolic_oh_pattern = Chem.MolFromSmarts("[c]O[H]")
      if mol.HasSubstructMatch(phenolic_oh_pattern):
        return True, "Contains cyclic terpene core and a phenolic hydroxyl."
    
    # 4. Check for fatty acid chains with specific functional groups (amides, esters) common in endocannabinoids
    fatty_acid_amide_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]=[CX3,CX4][CX3](=[OX1])[NX2]")
    if mol.HasSubstructMatch(fatty_acid_amide_pattern):
        return True, "Contains fatty acid chain with an amide."

    fatty_acid_ester_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]=[CX3,CX4][CX3](=[OX1])[OX2][CX4]")
    if mol.HasSubstructMatch(fatty_acid_ester_pattern):
         return True, "Contains fatty acid chain with an ester."

    # Negative constraint: if it is a peptide.
    peptide_pattern = Chem.MolFromSmarts("[NX2][CX3](=[OX1])[CX4]")
    if mol.HasSubstructMatch(peptide_pattern):
        return False, "Contains peptide bonds, unlikely to be a cannabinoid."


    return False, "Does not fit the criteria for a cannabinoid structure."