"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
"""
Classifies: unsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    An unsaturated fatty acyl-CoA results from the formal condensation of the thiol group
    of coenzyme A with the carboxy group of any unsaturated fatty acid.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is an unsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define Coenzyme A (CoA) substructure using SMARTS
    # Simplified pattern focusing on key functional groups
    coA_smarts = """
    NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](N2C=NC3=NC=NC=C23)[C@H](O)[C@@H]1OP(=O)(O)O
    """
    coA_mol = Chem.MolFromSmarts(coA_smarts)
    if coA_mol is None:
        return False, "Error parsing CoA SMARTS pattern"
    
    if not mol.HasSubstructMatch(coA_mol):
        return False, "No Coenzyme A moiety found"
    
    # Define thioester linkage pattern (C(=O)-S)
    thioester_smarts = "C(=O)S"
    thioester_mol = Chem.MolFromSmarts(thioester_smarts)
    if thioester_mol is None:
        return False, "Error parsing thioester SMARTS pattern"

    thioester_matches = mol.GetSubstructMatches(thioester_mol)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # For each thioester linkage, check if it's connected to CoA
    for match in thioester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon index
        sulfur_idx = match[1]      # Sulfur atom index

        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        # Check if sulfur is connected to CoA
        # Create a copy of the molecule to mark atoms
        mol_with_query = Chem.AddHs(mol)
        atom_props = {}
        for atom in mol_with_query.GetAtoms():
            atom_props[atom.GetIdx()] = False
        # Mark CoA atoms
        coa_matches = mol.GetSubstructMatch(coA_mol)
        for idx in coa_matches:
            atom_props[idx] = True
        # Check if sulfur is connected to CoA
        is_sulfur_connected_to_coa = atom_props[sulfur_idx]
        if not is_sulfur_connected_to_coa:
            continue  # Check next thioester linkage

        # Extract the fatty acyl chain (atoms connected to carbonyl carbon)
        fatty_acyl_atoms = Chem.GetMolFragmentAtoms(mol, [carbonyl_c_idx])
        fatty_acyl_bonds = []
        # Use BFS to traverse the acyl chain
        visited = set()
        to_visit = [carbonyl_c_idx]
        while to_visit:
            current_idx = to_visit.pop()
            if current_idx in visited:
                continue
            visited.add(current_idx)
            atom = mol.GetAtomWithIdx(current_idx)
            if atom.GetAtomicNum() in [6, 1]:  # Carbon or hydrogen
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    bond = mol.GetBondBetweenAtoms(current_idx, neighbor_idx)
                    fatty_acyl_bonds.append(bond.GetIdx())
                    to_visit.append(neighbor_idx)
        # Create a sub-molecule of the acyl chain
        acyl_chain = Chem.PathToSubmol(mol, fatty_acyl_bonds)
        # Check for unsaturation (C=C bonds)
        has_double_bond = False
        for bond in acyl_chain.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()
                if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                    has_double_bond = True
                    break
        if has_double_bond:
            return True, "Contains CoA moiety with unsaturated fatty acyl chain linked via thioester bond"
        else:
            return False, "Fatty acyl chain is saturated (no carbon-carbon double bonds)"

    return False, "No thioester linkage connected to CoA was found"