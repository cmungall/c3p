"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:57373 long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A long-chain fatty acyl-CoA results from the formal condensation of the thiol group
    of coenzyme A with the carboxy group of any long-chain (C13 to C22) fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for thioester group C(=O)S
    thioester_pattern = Chem.MolFromSmarts('C(=O)S')
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found"

    # Check for minimal CoA pattern (pantetheine moiety connected to S)
    coa_pattern = Chem.MolFromSmarts('SCCNC(=O)CCNC(=O)')
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "No CoA moiety found"

    # Now, for each thioester group, check if it's connected to CoA
    for match in thioester_matches:
        c_idx = match[0]  # Carbonyl carbon index
        o_idx = match[1]  # Carbonyl oxygen index
        s_idx = match[2]  # Sulfur index

        # Check if sulfur is part of the CoA moiety
        coa_connected = False
        for coa_match in coa_matches:
            if s_idx in coa_match:
                coa_connected = True
                break
        if not coa_connected:
            continue  # Sulfur atom not in CoA moiety

        # Get the acyl chain starting from the carbonyl carbon
        carbonyl_carbon = mol.GetAtomWithIdx(c_idx)
        # Exclude bonds to sulfur and oxygen
        acyl_neighbors = [nbr for nbr in carbonyl_carbon.GetNeighbors() if nbr.GetIdx() not in (o_idx, s_idx)]
        if not acyl_neighbors:
            continue  # No acyl chain attached to carbonyl carbon
        acyl_start_atom = acyl_neighbors[0]

        # Traverse the acyl chain and count the number of carbons
        visited = set()

        def count_acyl_carbons(atom, parent_idx):
            """
            Recursively counts the number of carbons in the acyl chain.
            """
            atom_idx = atom.GetIdx()
            if atom_idx in visited:
                return 0
            visited.add(atom_idx)
            count = 1 if atom.GetAtomicNum() == 6 else 0  # Count current atom if it's carbon
            for neighbor in atom.GetNeighbors():
                nbr_idx = neighbor.GetIdx()
                if nbr_idx != parent_idx and neighbor.GetAtomicNum() == 6:
                    bond = mol.GetBondBetweenAtoms(atom_idx, nbr_idx)
                    if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        count += count_acyl_carbons(neighbor, atom_idx)
            return count

        carbon_count = count_acyl_carbons(acyl_start_atom, carbonyl_carbon.GetIdx())

        # Exclude the carbonyl carbon from the count
        if 13 <= carbon_count <= 22:
            return True, f"Contains CoA moiety linked via thioester to a fatty acyl chain of length {carbon_count} carbons"
        else:
            return False, f"Acyl chain length is {carbon_count}, not in the range 13-22 carbons"

    return False, "Thioester group not connected to CoA moiety or acyl chain length not in range"