"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    A monoacyl-sn-glycerol 3-phosphate is a glycerol backbone with a single acyl chain 
    attached via an ester bond at position 1 or 2, and a phosphate group at position 3.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a monoacyl-sn-glycerol 3-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define phosphate attachment pattern: carbon connected via oxygen to phosphate group
    phosphate_attachment_pattern = Chem.MolFromSmarts("[C][O][P](=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_attachment_pattern)
    if not phosphate_matches:
        return False, "No glycerol phosphate backbone found"
    
    # Loop over phosphate matches
    for match in phosphate_matches:
        # match contains [C][O][P]
        c3_idx = match[0]
        o3_idx = match[1]
        p_idx = match[2]
        c3_atom = mol.GetAtomWithIdx(c3_idx)
        
        # Get neighbors of C3 besides O3
        c3_neighbors = [nbr.GetIdx() for nbr in c3_atom.GetNeighbors() if nbr.GetIdx() != o3_idx and nbr.GetAtomicNum() == 6]
        if len(c3_neighbors) < 1:
            continue  # Cannot find C2
        c2_idx = c3_neighbors[0]
        c2_atom = mol.GetAtomWithIdx(c2_idx)
        
        # Get neighbors of C2 besides C3
        c2_neighbors = [nbr.GetIdx() for nbr in c2_atom.GetNeighbors() if nbr.GetIdx() != c3_idx and nbr.GetAtomicNum() == 6]
        if len(c2_neighbors) < 1:
            continue  # Cannot find C1
        c1_idx = c2_neighbors[0]
        c1_atom = mol.GetAtomWithIdx(c1_idx)
        
        # Now check for ester linkage at C1 or C2
        # For each of C1 and C2, check if connected to an ester
        def has_ester_linkage(c_idx):
            atom = mol.GetAtomWithIdx(c_idx)
            for bond in atom.GetBonds():
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    o_atom = nbr
                    # Check if oxygen is connected to a carbonyl carbon (C=O)
                    for obond in o_atom.GetBonds():
                        o_nbr = obond.GetOtherAtom(o_atom)
                        if o_nbr.GetIdx() == c_idx:
                            continue
                        if o_nbr.GetAtomicNum() == 6 and obond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            # Check for double bond (C=O)
                            carbonyl_carbon = o_nbr
                            has_c_double_o = False
                            for cbond in carbonyl_carbon.GetBonds():
                                if cbond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                    o_double = cbond.GetOtherAtom(carbonyl_carbon)
                                    if o_double.GetAtomicNum() == 8:
                                        has_c_double_o = True
                                        break
                            if has_c_double_o:
                                return True
            return False
        
        def has_hydroxyl(c_idx):
            atom = mol.GetAtomWithIdx(c_idx)
            for bond in atom.GetBonds():
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    if len(nbr.GetNeighbors()) == 1:
                        return True
            return False
        
        ester_positions = []
        hydroxyl_positions = []
        for c_idx in [c1_idx, c2_idx]:
            if has_ester_linkage(c_idx):
                ester_positions.append(c_idx)
            elif has_hydroxyl(c_idx):
                hydroxyl_positions.append(c_idx)
        
        if len(ester_positions) == 1 and len(hydroxyl_positions) == 1:
            return True, "Contains glycerol backbone with phosphate group at position 3, and single acyl chain attached via ester linkage at position 1 or 2"
    
    return False, "Does not match monoacyl-sn-glycerol 3-phosphate structure"