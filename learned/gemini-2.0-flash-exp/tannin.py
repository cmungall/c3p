"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are astringent polyphenolic vegetable principles, chiefly complex glucosides of catechol and pyrogallol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    score = 0 # Start with score of zero
    reason = ""

    # 1. Check for linked aromatic rings
    linked_benzene_pattern = Chem.MolFromSmarts("c1ccccc1(-c2ccccc2)")
    linked_benzene_matches = mol.GetSubstructMatches(linked_benzene_pattern)
    if len(linked_benzene_matches) > 0:
        score += len(linked_benzene_matches) * 2  # Award points for linked aromatic rings.

    #Check for multiple aromatic rings - basic check
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    if len(benzene_matches) >= 2:
        score += 2
    elif len(benzene_matches) > 0: #Only one aromatic ring is not enough, but we give it a point
        score += 1

    # 2. Check for catechol or pyrogallol units and variations
    catechol_pattern_1 = Chem.MolFromSmarts("c1c(O)c(O)cc1")
    catechol_pattern_2 = Chem.MolFromSmarts("c1c(O)c(O)c([OX2])c1") #Ether variations
    catechol_matches = mol.GetSubstructMatches(catechol_pattern_1) + mol.GetSubstructMatches(catechol_pattern_2)

    pyrogallol_pattern_1 = Chem.MolFromSmarts("c1c(O)c(O)c(O)cc1")
    pyrogallol_pattern_2 = Chem.MolFromSmarts("c1c(O)c(O)c(O)c([OX2])c1")
    pyrogallol_matches = mol.GetSubstructMatches(pyrogallol_pattern_1) + mol.GetSubstructMatches(pyrogallol_pattern_2)

    if len(catechol_matches) > 0:
        score += len(catechol_matches) * 2
    if len(pyrogallol_matches) > 0:
        score += len(pyrogallol_matches) * 3
    
    # 3. Check for glycosidic bonds (C-O-C with an adjacent 5 or 6 membered ring containing O)
    glycosidic_pattern_6 = Chem.MolFromSmarts("[CX4]-[OX2]-[C]1[CH2][CH2][CH][CH][CH][O]1") #6 membered ring
    glycosidic_pattern_5 = Chem.MolFromSmarts("[CX4]-[OX2]-[C]1[CH2][CH2][CH][CH][O]1") #5 membered ring
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern_6) + mol.GetSubstructMatches(glycosidic_pattern_5)
    score += len(glycosidic_matches)
    
    # 4. Check for catechin core (simplified)
    catechin_core_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c(C[C@H]2[C@H](O)[CH](O)c3c(O)cc(O)cc32)c1")
    catechin_core_matches = mol.GetSubstructMatches(catechin_core_pattern)
    if len(catechin_core_matches) > 0:
        score += len(catechin_core_matches) * 4 #High score for a catechin core.

    # 5. Check for Gallic acid component
    gallic_acid_pattern = Chem.MolFromSmarts("c1c(O)c(O)c(C(=O)O)c(O)c1")
    gallic_matches = mol.GetSubstructMatches(gallic_acid_pattern)
    score += len(gallic_matches)


    # 6. Check molecular weight and size
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if 500 <= mol_wt <= 2000 : #typical range for a tannin
        score += 2
    elif mol_wt > 2000: #Penalise too large molecules
        score -= 1
    
    num_atoms = mol.GetNumAtoms()
    if num_atoms >= 30: #Rough check for size
      score += 1


    # Final classification based on score.
    if score >= 8:
        return True, "Molecule contains multiple linked aromatic rings, hydroxyl groups, and catechol/pyrogallol units, with catechin core, glycosidic bonds and appropriate MW"
    elif score >= 5:
        return False, "Molecule has some tannin-like features, but not enough"
    else:
         return False, "Molecule does not resemble a tannin"