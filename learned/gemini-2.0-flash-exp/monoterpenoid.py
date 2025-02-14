"""
Classifies: CHEBI:25409 monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
         bool: True if molecule is a monoterpenoid, False otherwise
         str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Hints for molecular weight and carbon count, less strict now
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Isoprene units check
    isoprene_pattern = Chem.MolFromSmarts("CC(=C)C") # Basic isoprene unit
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    
    
    # More specific substructure patterns
    menthane_pattern = Chem.MolFromSmarts("[C]1([C])[C]([C])([C])[C]([C])[C]1")
    pinane_pattern = Chem.MolFromSmarts("C1[C]2[C]([C]1([C])C)[C]([C])([C])CC2")
    bornane_pattern = Chem.MolFromSmarts("C1[C]2[C]([C]1([C])C)[C]([C])([C])CC2")
    thujane_pattern = Chem.MolFromSmarts("C1[C]2[C]([C]1([C])C)[C]([C])([C])CC2")

    p_menthene_pattern1 = Chem.MolFromSmarts("C[C]1([C])[C]([C])([C])[C]([C])=[C]1")
    p_menthene_pattern2 = Chem.MolFromSmarts("C[C]1([C])=[C]([C])[C]([C])[C]([C])[C]1")
    
    # Include more complex variations
    menthene_pattern_carbonyl = Chem.MolFromSmarts("[C]1([C])[C]([C])([C])[C](=[O])[C]1")
    menthene_pattern_hydroxyl = Chem.MolFromSmarts("[C]1([C])[C]([C])([C])[C]([OH])[C]1")    
    
    bicyclic_pinene = Chem.MolFromSmarts("C1[C]2[C](C1([C])C)[C]([C])=[C]C2")
    bicyclic_thujene = Chem.MolFromSmarts("C1[C]2[C](C1([C])C)[C]([C])=[C]C2")
    bicyclic_bornane = Chem.MolFromSmarts("C1[C]2[C](C1([C])C)[C]([C])([C])C2")
    
    
    if (mol.HasSubstructMatch(menthane_pattern) or
            mol.HasSubstructMatch(pinane_pattern) or
            mol.HasSubstructMatch(bornane_pattern) or
            mol.HasSubstructMatch(thujane_pattern) or
            mol.HasSubstructMatch(p_menthene_pattern1) or
            mol.HasSubstructMatch(p_menthene_pattern2) or
            mol.HasSubstructMatch(menthene_pattern_carbonyl) or
            mol.HasSubstructMatch(menthene_pattern_hydroxyl) or
            mol.HasSubstructMatch(bicyclic_pinene) or
            mol.HasSubstructMatch(bicyclic_thujene) or
            mol.HasSubstructMatch(bicyclic_bornane)
            ):
            
        if len(isoprene_matches) >= 1:
                return True, "Matches monoterpenoid substructure and contains one or more isoprene unit(s)"
        else:
             return False, "Matches monoterpenoid substructure, but is missing isoprene unit(s)"
            
    elif len(isoprene_matches) >= 2 and (c_count >= 8 and c_count <= 15):
        return True, "Contains two or more isoprene units, and carbon count suggest it could be a monoterpenoid"

    # Ring Check: allow for fused bicyclics (relax constraints)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 0 :
        found_6_member = False
        found_5_member = False
        for ring in ring_info.AtomRings():
            if len(ring) == 6:
              found_6_member = True
            if len(ring) == 5:
                found_5_member = True
        if not (found_6_member or found_5_member):
           if  not (mol.HasSubstructMatch(bicyclic_pinene) or mol.HasSubstructMatch(bicyclic_thujene) or mol.HasSubstructMatch(bicyclic_bornane)):
               return False, "Monoterpenoids usually contain a 5 or 6-membered ring or a bicyclic fused ring"

    if  (c_count >= 8 and c_count <= 15) and (mol_wt >=120 and mol_wt <= 400):
            if  len(isoprene_matches) >= 1:
               return True, "Carbon count, molecular weight, and isoprene matches suggest it could be a monoterpenoid"

    return False, "Does not match common monoterpenoid substructures or characteristics"