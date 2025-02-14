"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA is an unsaturated fatty acyl-CoA with a double bond between positions 2 and 3
    of the acyl chain attached via a thioester linkage to Coenzyme A.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for Coenzyme A (CoA) moiety
    coa_smarts = '[OH][C@H](COP(O)(=O)OP(O)(=O)OCC1OC(CO[P](O)(=O)O)C(O)C(O)C1O)n1cnc2c(N)ncnc12'
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A (CoA) moiety not found"
    
    # Define SMARTS pattern for thioester linkage: S-C(=O)-C
    thioester_smarts = 'SC(=O)C'
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"
    
    # Assuming only one acyl chain attached via thioester linkage
    thioester_match = thioester_matches[0]
    
    # Get the carbon atom of the acyl chain (first carbon after the carbonyl)
    acyl_c1_idx = thioester_match[2]  # Index of the first carbon in the acyl chain
    acyl_chain_atoms = [acyl_c1_idx]
    acyl_chain_bonds = []
    
    # Traverse the acyl chain to collect atoms and bonds
    mol = Chem.AddHs(mol)
    atom_queue = [acyl_c1_idx]
    visited_atoms = set()
    while atom_queue:
        atom_idx = atom_queue.pop(0)
        if atom_idx in visited_atoms:
            continue
        visited_atoms.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() != 'C':
            continue  # Stop traversal if not carbon
        for bond in atom.GetBonds():
            nbr_atom = bond.GetOtherAtom(atom)
            nbr_idx = nbr_atom.GetIdx()
            if nbr_idx in visited_atoms:
                continue
            if nbr_atom.GetAtomicNum() == 6:
                acyl_chain_atoms.append(nbr_idx)
                acyl_chain_bonds.append(bond)
                atom_queue.append(nbr_idx)
            elif nbr_atom.GetAtomicNum() == 1:
                continue  # Ignore hydrogens
            else:
                continue  # Stop at heteroatoms
    
    # Sort atoms by their indices to represent the chain order
    acyl_chain_atoms.sort()
    
    # Check for double bond between C2 and C3 (positions start from acyl_c1_idx)
    double_bond_found = False
    for bond in acyl_chain_bonds:
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        bond_order = bond.GetBondTypeAsDouble()
        if bond_order == 2.0:
            # Adjust indices to start from 1
            pos1 = acyl_chain_atoms.index(begin_idx) + 1
            pos2 = acyl_chain_atoms.index(end_idx) + 1
            if {pos1, pos2} == {2, 3}:
                double_bond_found = True
                break
    
    if not double_bond_found:
        return False, "No double bond between positions 2 and 3 in the acyl chain"
    
    return True, "Molecule is a 2-enoyl-CoA with a double bond between positions 2 and 3"