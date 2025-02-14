"""
Classifies: CHEBI:64611 ether lipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid has at least one fatty acid chain attached to the glycerol backbone via an ether linkage (C-O-C).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached, and no carbonyl)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4;!$(C=O)][CHX4;!$(C=O)][CH2X4;!$(C=O)]")
    if not mol.HasSubstructMatch(glycerol_pattern):
         glycerol_pattern = Chem.MolFromSmarts("[CH2X4;!$(C=O)]=[CHX3;!$(C=O)][CH2X4;!$(C=O)]")
         if not mol.HasSubstructMatch(glycerol_pattern):
            return False, "No glycerol backbone found"

    # Get the indices of the glycerol carbons
    glycerol_carbons = []
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if glycerol_matches:
         glycerol_carbons = glycerol_matches[0]
    else:
         return False, "Glycerol substructure not found"

    # Look for at least one ether group connected to a glycerol carbon (C-O-C),
    # with the glycerol carbon connected to the ether oxygen, and also to a long chain
    ether_pattern = Chem.MolFromSmarts("[CH2X4,CHX4;!$(C=O)]-[OX2]-[CX4]-[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    ether_pattern2 = Chem.MolFromSmarts("[CX3]=[CX3]-[OX2]-[CX4]-[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")

    ether_matches = mol.GetSubstructMatches(ether_pattern)
    ether_matches2 = mol.GetSubstructMatches(ether_pattern2)


    has_valid_ether = False
    if ether_matches:
      for match in ether_matches:
        #check if the oxygen is connected to a glycerol carbon
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'O':
                neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
                glycerol_neighbors = [n for n in neighbors if n in glycerol_carbons]
                if glycerol_neighbors:
                   has_valid_ether = True
                   break
        if has_valid_ether:
            break
    if not has_valid_ether and ether_matches2:
       for match in ether_matches2:
        #check if the oxygen is connected to a glycerol carbon
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'O':
                neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
                glycerol_neighbors = [n for n in neighbors if n in glycerol_carbons]
                if glycerol_neighbors:
                   has_valid_ether = True
                   break
        if has_valid_ether:
            break

    if not has_valid_ether:
        return False, "No valid ether linkage found on glycerol"

    # Check molecular weight - lipids typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for ether lipid"


    return True, "Contains glycerol backbone with at least one fatty acid chain attached via ether bond"