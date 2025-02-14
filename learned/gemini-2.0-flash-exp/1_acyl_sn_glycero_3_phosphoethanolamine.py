"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors


def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    These molecules have a glycerol backbone esterified at the sn-1 position with a fatty acid,
    and at the sn-3 position with a phosphoethanolamine group and (R) configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the glycerol backbone with the correct stereochemistry using SMARTS
    # [C@H](O)(CO) is the correct stereochemistry.
    glycerol_pattern = Chem.MolFromSmarts("[C@H](O)(CO)")
    glycerol_match = mol.GetSubstructMatch(glycerol_pattern)

    if not glycerol_match:
        return False, "No glycerol backbone with correct stereochemistry found"

    # Check for the presence of a fatty acid ester group at sn-1
    # This is one of the carbon attached to the glycerol backbone
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not ester_matches:
         return False, f"No ester group at sn-1 glycerol position"
    
    # Check for phosphoethanolamine group
    phosphoethanolamine_pattern = Chem.MolFromSmarts("COP(=O)(O)OCCN")
    phosphoethanolamine_match = mol.GetSubstructMatch(phosphoethanolamine_pattern)

    if not phosphoethanolamine_match:
        return False, "No phosphoethanolamine group found"

    # Verify position of the phosphoethanolamine group on the glycerol
    # The phosphorus atom should be connected to the CO group of the glycerol
    # which is connected to the [C@H] carbon
    
    glycerol_carbon_indices = glycerol_match
    ester_carbon_indices = []
    for match in ester_matches:
        ester_carbon_indices.append(match[2])
    
    if len(glycerol_carbon_indices) != 3:
        return False, "Incorrect number of atoms in the glycerol backbone"
    
    # The glycerol carbon connected to the ester should not be connected to the phosphorus
    
    glycerol_carbon_ester = -1
    for index in ester_carbon_indices:
        if mol.GetAtomWithIdx(index).GetBonds()[0].GetOtherAtomIdx(index) == glycerol_carbon_indices[0]:
            glycerol_carbon_ester = index
    
    if glycerol_carbon_ester == -1:
         return False, "Ester group not correctly attached to the glycerol backbone"

    
    phosphorus_index = -1
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:
           phosphorus_index = atom.GetIdx()
           break

    if phosphorus_index == -1:
        return False, "No phosphorus atom found"

    
    phosphoethanolamine_attached_carbon = -1
    for atom in mol.GetAtomWithIdx(phosphorus_index).GetNeighbors():
        if atom.GetAtomicNum() == 8:
           for atom2 in mol.GetAtomWithIdx(atom.GetIdx()).GetNeighbors():
                if atom2.GetAtomicNum() == 6:
                    phosphoethanolamine_attached_carbon = atom2.GetIdx()
                    break
           if phosphoethanolamine_attached_carbon != -1:
             break

    if phosphoethanolamine_attached_carbon == -1:
        return False, "Phosphate group not correctly attached"
        
    if not mol.GetAtomWithIdx(phosphoethanolamine_attached_carbon).IsInRing():
        if mol.GetAtomWithIdx(phosphoethanolamine_attached_carbon).GetBonds()[0].GetOtherAtomIdx(phosphoethanolamine_attached_carbon) == glycerol_carbon_indices[1] :
             pass
        else:
            return False, "Phosphate group not attached to the correct carbon"
    else:
        return False, "Phosphate group must be attached to the exocyclic carbon of glycerol"
    # Check for a fatty acid chain (long carbon chain)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
         return False, "Missing fatty acid chain"
        
    #Check for C(=O)
    
    carbonyl_pattern = Chem.MolFromSmarts("C=O")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    
    if len(carbonyl_matches) == 0:
        return False, "No carbonyl groups found"
    
    # Verify if the carbonyl carbon is attached to the glycerol carbon
    carbonyl_attached = False
    for match in carbonyl_matches:
        carbonyl_carbon = match[0]
        for neighbor in mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors():
            if neighbor.GetIdx() == glycerol_carbon_ester:
                carbonyl_attached = True
                break

    if not carbonyl_attached:
       return False, "Carbonyl group not attached to the glycerol backbone"

    # Check molecular weight - typically > 400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
         return False, "Molecular weight too low for a phosphoethanolamine"
    
    #Check for exactly one glycerol, ester, and phosphoethanolamine group
    if len(glycerol_match) == 0:
        return False, "Incorrect number of glycerol groups"
    if len(ester_matches) == 0:
        return False, "Incorrect number of ester groups"
    if len(phosphoethanolamine_match) == 0:
        return False, "Incorrect number of phosphoethanolamine groups"


    return True, "Matches 1-acyl-sn-glycero-3-phosphoethanolamine criteria"