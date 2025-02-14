"""
Classifies: CHEBI:33447 phospho sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is a monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Monosaccharide ring patterns (generic 5 or 6 membered ring with O and C)
    # This allows for substitutions on the ring
    monosaccharide_pattern = Chem.MolFromSmarts("[CX4,CX3][CX4,CX3][CX4,CX3][OX2][CX4,CX3][CX4,CX3]")
    monosaccharide_pattern2 = Chem.MolFromSmarts("[CX4,CX3][CX4,CX3][CX4,CX3][CX4,CX3][OX2][CX4,CX3]")

    if not mol.HasSubstructMatch(monosaccharide_pattern) and not mol.HasSubstructMatch(monosaccharide_pattern2):
        return False, "No monosaccharide ring found"
    
    # Phosphate group pattern (ester)
    phosphate_pattern = Chem.MolFromSmarts("[OX2][P](=[OX1])([OX2])([OX2])")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"
    
    # Check for direct linkage of phosphate to a monosaccharide oxygen
    linked = False
    for match in phosphate_matches:
        for p_idx in match: #iterate through atoms in the phosphate substructure
            p_atom = mol.GetAtomWithIdx(p_idx)
            if p_atom.GetAtomicNum() == 8: # get only the oxygens from the phosphate
                for neighbor in p_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6: # if bonded to a carbon
                        for n2 in neighbor.GetNeighbors():
                           if n2.GetAtomicNum() == 8 and n2.IsInRing(): #Check if the oxygen bonded to that carbon is in a ring
                            linked = True
                            break
                           elif n2.GetAtomicNum() == 8 and not n2.IsInRing() and len(n2.GetNeighbors()) == 1: #Check if the oxygen bonded to that carbon is an alcohol
                            linked = True
                            break
                    if linked:
                       break
                if linked:
                  break
        if linked:
            break
    
    if not linked:
        return False, "Phosphate not directly linked to monosaccharide"

    return True, "Contains a monosaccharide with a phosphate group attached via ester bond"