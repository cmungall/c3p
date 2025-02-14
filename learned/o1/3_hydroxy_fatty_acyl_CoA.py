"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA is a hydroxy fatty acyl-CoA that results from the formal condensation 
    of the thiol group of coenzyme A with the carboxy group of any 3-hydroxy fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for thioester linkage C(=O)S
    thioester_pattern = Chem.MolFromSmarts('C(=O)S')
    matches = mol.GetSubstructMatches(thioester_pattern)

    if not matches:
        return False, "No thioester linkage found"

    # Iterate over all thioester matches
    for match in matches:
        carbonyl_carbon_idx = match[0]
        sulfur_idx = match[2]

        carbonyl_carbon = mol.GetAtomWithIdx(carbonyl_carbon_idx)
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)

        # Check if sulfur is connected to CoA moiety (simplified check)
        # We assume that sulfur connected to nitrogen atoms is indicative of CoA
        sulfur_neighbors = sulfur_atom.GetNeighbors()
        is_connected_to_coa = False
        for neighbor in sulfur_neighbors:
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                # Check if this carbon is connected to nitrogen atoms
                carbon_neighbors = neighbor.GetNeighbors()
                nitrogen_count = sum(1 for atom in carbon_neighbors if atom.GetAtomicNum() == 7)
                if nitrogen_count >= 1:
                    is_connected_to_coa = True
                    break
        if not is_connected_to_coa:
            continue  # Not connected to CoA, check next match

        # Identify C2 atom in acyl chain (attached to carbonyl carbon and not sulfur)
        c2_atom = None
        for neighbor in carbonyl_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != sulfur_idx:
                c2_atom = neighbor
                break
        if c2_atom is None:
            continue  # No acyl chain found

        # Identify C3 atom (next carbon in the chain)
        c3_atom = None
        for neighbor in c2_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != carbonyl_carbon_idx:
                c3_atom = neighbor
                break
        if c3_atom is None:
            continue  # Acyl chain too short

        # Check for hydroxy group on C3 atom
        has_hydroxy = False
        for neighbor in c3_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen atom
                if neighbor.GetDegree() == 1:
                    # Check bond type (should be single bond)
                    bond = mol.GetBondBetweenAtoms(c3_atom.GetIdx(), neighbor.GetIdx())
                    if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        has_hydroxy = True
                        break
        if has_hydroxy:
            return True, "Contains 3-hydroxy group on acyl chain attached to CoA via thioester linkage"

    return False, "Does not contain 3-hydroxy group on acyl chain attached to CoA via thioester linkage"