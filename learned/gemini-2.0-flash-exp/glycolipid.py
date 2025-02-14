"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is a lipid (fatty acids, diacylglycerol, or ceramide) linked to a carbohydrate (mono-, di-, or trisaccharide) through a glycosidic bond.
    Some glycolipids don't have glycerol or sphingosine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """

    # 1. Basic Checks
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Sugar Ring and Glycosidic Linkage
    sugar_ring_pattern = Chem.MolFromSmarts("[OX2]([CX4])[CX4][CX4][CX4][CX4]") # A 6-membered ring
    sugar_ring_pattern2 = Chem.MolFromSmarts("[OX2]([CX4])[CX4][CX4][CX4]") # A 5-membered ring
    sugar_matches = mol.GetSubstructMatches(sugar_ring_pattern)
    sugar_matches2 = mol.GetSubstructMatches(sugar_ring_pattern2)
    if not sugar_matches and not sugar_matches2:
        return False, "No sugar ring found"

    glycosidic_pattern = Chem.MolFromSmarts("[OX2]([CX4])-[CX4]")
    glycosidic_matches = []
    for match in mol.GetSubstructMatches(glycosidic_pattern):
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 8:
                # Find if one of the carbon connected is part of the sugar ring
                connected_carbons = []
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:
                       connected_carbons.append(neighbor.GetIdx())
                
                is_sugar_carbon = False
                for sugar_match in sugar_matches + sugar_matches2:
                    if connected_carbons[0] in sugar_match or connected_carbons[1] in sugar_match:
                         is_sugar_carbon = True
                         break
                if is_sugar_carbon:
                  glycosidic_matches.append(match)


    if not glycosidic_matches:
        return False, "No glycosidic linkage from sugar ring found"

    # 3. Lipid Component Detection (Fatty acids, Glycerol, Ceramides)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)

    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern)

    ceramide_pattern = Chem.MolFromSmarts("[CX4][NX3][CX4]~[CX4]")
    has_ceramide = mol.HasSubstructMatch(ceramide_pattern)
    
    if not fatty_acid_matches and not has_glycerol and not has_ceramide:
        return False, "No lipid component detected"

    # Check for rotatable bonds for fatty chains
    if fatty_acid_matches:
        n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
        if n_rotatable < 5:
              return False, "Chains too short to be fatty acids"

    # 4. Molecular Weight Check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for glycolipid"
    
    # 5. Atom count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:
        return False, "Too few carbons for glycolipid"
    if o_count < 3:
         return False, "Too few oxygens for glycolipid"

    return True, "Glycosidic linkage, carbohydrate and lipid chains are present."