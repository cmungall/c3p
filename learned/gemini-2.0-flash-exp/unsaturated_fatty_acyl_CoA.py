"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
"""
Classifies: unsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    An unsaturated fatty acyl-CoA is a molecule containing a coenzyme A moiety linked to an
    unsaturated fatty acid via a thioester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for charges
    if any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms()):
        return False, "Molecule contains charges"
    
    # Define the CoA substructure using SMARTS
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Check for the thioester bond (-C(=O)-S-)
    thioester_pattern = Chem.MolFromSmarts("C(=O)SC")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) != 1:
        return False, f"Found {len(thioester_matches)} thioester groups, need exactly 1"

    # Find the carbon atom of the carbonyl group in thioester
    thioester_c_index = thioester_matches[0][0]
    
    # Find the fatty acid chain connected to the thioester
    fatty_acid_atoms = set()
    queue = [thioester_c_index]
    visited = {thioester_c_index}

    while queue:
        current_idx = queue.pop(0)
        current_atom = mol.GetAtomWithIdx(current_idx) # get atom object
        for neighbor in current_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited:
                if neighbor.GetSymbol() == 'C':
                    fatty_acid_atoms.add(neighbor_idx)
                    queue.append(neighbor_idx)
                visited.add(neighbor_idx)
            
    # Check for at least one double bond (C=C) in the acyl chain (exclude double bonds in the CoA)
    double_bond_pattern = Chem.MolFromSmarts("[CH]=[CH]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    acyl_double_bonds = 0

    for match in double_bond_matches:
        atom1_idx = mol.GetAtomWithIdx(match[0]).GetIdx()
        atom2_idx = mol.GetAtomWithIdx(match[1]).GetIdx()
        if atom1_idx in fatty_acid_atoms and atom2_idx in fatty_acid_atoms:
            acyl_double_bonds += 1

    if acyl_double_bonds < 1:
      return False, "No double bond found in the fatty acid chain"
    
    # Check the number of carbons in the fatty acid chain.
    chain_carbon_count = len(fatty_acid_atoms)
    if chain_carbon_count < 4 :
      return False, f"Fatty acid chain too short, found {chain_carbon_count} carbons"
      
    # Check molecular weight - Acyl-CoAs typically > 600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for unsaturated fatty acyl-CoA"

    return True, "Contains CoA moiety linked to an unsaturated fatty acid via a thioester bond"