"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    A monoacyl-sn-glycerol 3-phosphate has a glycerol backbone, a phosphate group at the 3-position,
    and a single fatty acid chain (acyl group) at either position 1 or 2 of the glycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacyl-sn-glycerol 3-phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sn-glycerol-3-phosphate core with correct stereochemistry
    # Using multiple SMARTS to cover different representations
    glycerol_phosphate_patterns = [
        Chem.MolFromSmarts("[CH2][C@H](O)COP(=O)(O)(O)"), # C1-C2(R)-C3-O-P
        Chem.MolFromSmarts("[CH2][C@H](O)C[O]P(=O)(O)(O)"), # C1-C2(R)-C3-[O]-P
        Chem.MolFromSmarts("[CH2][C@H](O)C[OP](=O)(O)O"),  #C1-C2(R)-C3-O-P
        Chem.MolFromSmarts("[CH2][C@H](O)C[O][P](O)(O)=O"), # C1-C2(R)-C3-[O]-P
        Chem.MolFromSmarts("[CH2][C@H](O)COP([O-])([O-])=O"), # C1-C2(R)-C3-O-P (charged form)
        Chem.MolFromSmarts("[CH2][C@H](O)C[O][P]([O-])([O-])=O"),  #C1-C2(R)-C3-[O]-P (charged form)
        Chem.MolFromSmarts("[CH2]C[C@@H](O)OP(=O)(O)(O)"), # C3-C2(S)-C1-O-P (enantiomer - not sn)
        Chem.MolFromSmarts("C[C@@H](O)C[O]P(=O)(O)(O)") # C2(S)-C3-[O]-P (enantiomer - not sn)
    ]
    
    glycerol_phosphate_matches = []
    for pattern in glycerol_phosphate_patterns:
        matches = mol.GetSubstructMatches(pattern)
        if matches:
           glycerol_phosphate_matches = matches
           break # Stop on the first match

    if not glycerol_phosphate_matches:
        return False, "No sn-glycerol-3-phosphate core found"

    # Verify the chirality of the central carbon is S (CCW)
    glycerol_carbon_index = glycerol_phosphate_matches[0][1]
    chiral_center = mol.GetAtomWithIdx(glycerol_carbon_index)
    if chiral_center.GetChiralTag() != Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
        return False, "Central carbon is not S chiral"

    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    #check ester connected to either position 1 or 2. Not the phosphate
    ester_match = ester_matches[0]
    carbonyl_carbon_index = ester_match[1]
    
    glycerol_phosphate_match = glycerol_phosphate_matches[0]
    
    glycerol_carbon1_index = glycerol_phosphate_match[0] #C1
    glycerol_carbon2_index = glycerol_phosphate_match[1] #C2
    
    ester_connected = False
    
    #Get the neighboring oxygen for each glycerol carbon
    glycerol_carbon1_oxygen_index = None
    for neighbor in mol.GetAtomWithIdx(glycerol_carbon1_index).GetNeighbors():
        if neighbor.GetAtomicNum() == 8:
           glycerol_carbon1_oxygen_index = neighbor.GetIdx()
           break
    
    glycerol_carbon2_oxygen_index = None
    for neighbor in mol.GetAtomWithIdx(glycerol_carbon2_index).GetNeighbors():
        if neighbor.GetAtomicNum() == 8:
            glycerol_carbon2_oxygen_index = neighbor.GetIdx()
            break

    if glycerol_carbon1_oxygen_index is not None:
        for neighbor in mol.GetAtomWithIdx(glycerol_carbon1_oxygen_index).GetNeighbors():
            if neighbor.GetIdx() == carbonyl_carbon_index:
                ester_connected = True
                break
        
    if not ester_connected and glycerol_carbon2_oxygen_index is not None:
      for neighbor in mol.GetAtomWithIdx(glycerol_carbon2_oxygen_index).GetNeighbors():
        if neighbor.GetIdx() == carbonyl_carbon_index:
            ester_connected = True
            break
           
    if not ester_connected:
       return False, "Ester not connected to position 1 or 2 of glycerol"
       
    # Check for long carbon chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
      return False, "Missing acyl chain"
        
    return True, "Contains a glycerol backbone, a phosphate group at position 3, and a single acyl group at position 1 or 2"