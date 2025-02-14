"""
Classifies: CHEBI:76579 triradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.
    A triradylglycerol has a glycerol backbone with three radyl groups attached.
    Radyl groups can be acyl, alkyl or alk-1-enyl.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triradylglycerol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Glycerol backbone check: same as for triglycerides.
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # 2. Check for three radyl groups - separate checks for acyl, alkyl, and alk-1-enyl
    acyl_pattern = Chem.MolFromSmarts("[OX2]-[CX3](=[OX1])")
    alkyl_pattern = Chem.MolFromSmarts("[OX2]-[CX4]")
    alkenyl_pattern = Chem.MolFromSmarts("[OX2]-[CX3]=[CX3]")


    radyl_matches = []

    glycerol_atoms = mol.GetSubstructMatches(glycerol_pattern)[0]
    
    for atom in glycerol_atoms:
        attached_atoms = []
        
        glycerol_atom = mol.GetAtomWithIdx(atom)
        for neighbor in glycerol_atom.GetNeighbors():
            if neighbor.GetIdx() not in glycerol_atoms:
                attached_atoms.append(neighbor)

        for attached_atom in attached_atoms:
            is_radyl = False
            
            #check acyl
            acyl_match = mol.GetSubstructMatches(acyl_pattern)
            for match in acyl_match:
                if attached_atom.GetIdx() == match[0]:
                    is_radyl = True
                    break
            if is_radyl:
               radyl_matches.append(attached_atom.GetIdx())
               continue
            
            #check alkyl
            alkyl_match = mol.GetSubstructMatches(alkyl_pattern)
            for match in alkyl_match:
                if attached_atom.GetIdx() == match[0]:
                    is_radyl = True
                    break
            if is_radyl:
                radyl_matches.append(attached_atom.GetIdx())
                continue

            #check alkenyl
            alkenyl_match = mol.GetSubstructMatches(alkenyl_pattern)
            for match in alkenyl_match:
                if attached_atom.GetIdx() == match[0]:
                    is_radyl = True
                    break
            if is_radyl:
                radyl_matches.append(attached_atom.GetIdx())
                continue
                
            if not is_radyl:
               return False, "Not all positions have valid radyl groups"
                
    
    if len(radyl_matches) != 3:
       return False, f"Found {len(radyl_matches)} radyl groups, need exactly 3"

    # 3. Check for long carbon chains
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, f"Too few carbons ({c_count}) in chains"

    # 4. Check number of oxygens: should be 6
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 6:
        return False, "Must have exactly 6 oxygens"
    
    # 5. Check presence of other elements, apart from C, H and O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1,6,8]:
            return False, "Molecule contains elements other than C, H and O"

    return True, "Contains glycerol backbone with three radyl groups (acyl/alkyl/alk-1-enyl)"