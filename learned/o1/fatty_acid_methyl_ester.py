"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
from rdkit import Chem

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is the carboxylic ester obtained by the formal condensation 
    of a fatty acid with methanol.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the methyl ester functional group pattern
    methyl_ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(methyl_ester_pattern)
    if not ester_matches:
        return False, "No methyl ester functional group found"
    
    # Analyze each methyl ester group found
    for match in ester_matches:
        carbonyl_c_idx = match[0]
        carbonyl_o_idx = match[1]
        ester_o_idx = match[2]
        methyl_c_idx = match[3]
        
        # Ensure the methyl group is terminal
        methyl_c = mol.GetAtomWithIdx(methyl_c_idx)
        if methyl_c.GetDegree() != 1:
            continue  # Not a methyl group
        
        # Get the carbonyl carbon
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
        
        # Exclude the ester group atoms from the acyl chain
        excluded_atoms = set(match)
        
        # Traverse the acyl chain starting from the carbonyl carbon
        acyl_chain_atoms = set()
        atoms_to_visit = [nbr for nbr in carbonyl_c.GetNeighbors() if nbr.GetIdx() not in excluded_atoms]
        while atoms_to_visit:
            atom = atoms_to_visit.pop()
            idx = atom.GetIdx()
            if idx in acyl_chain_atoms or idx in excluded_atoms:
                continue  # Already visited or excluded
            acyl_chain_atoms.add(idx)
            # Allow typical elements: C, H, O, N, S
            if atom.GetAtomicNum() not in (6, 1, 8, 7, 16):
                return False, f"Element {atom.GetSymbol()} not typical in fatty acid chain"
            # Add neighbors to visit
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx not in acyl_chain_atoms and nbr_idx not in excluded_atoms:
                    atoms_to_visit.append(nbr)
        
        # Count the number of carbons in the acyl chain
        c_count = sum(1 for idx in acyl_chain_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if c_count < 4:
            return False, "Acyl chain too short to be a fatty acid"
        
        # Passed all checks
        return True, "Contains methyl ester group with appropriate fatty acid chain"
    
    return False, "No suitable fatty acid methyl ester pattern found"