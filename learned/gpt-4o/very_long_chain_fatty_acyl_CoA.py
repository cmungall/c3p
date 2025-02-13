"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA 
    based on its SMILES string. Very long-chain fatty acyl-CoAs 
    have fatty acyl groups with chain lengths greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA,
              False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES into RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify Coenzyme A part by locating substructures.
    # Assume CoA is a fixed known structure in the molecule with common residues.
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    coa_match = mol.GetSubstructMatch(coa_pattern)
    if not coa_match:
        return False, "No Coenzyme A (CoA) core structure found"
    
    # The position of thioester linkage, carbon for acyl chain, is directly connected to sulfur.
    thioester_carbon = None
    for atom_idx in coa_match:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() == 'S':
            # Iterate over neighbors to find the relevant carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    thioester_carbon = neighbor
                    break

    if thioester_carbon is None:
        return False, "Thioester bond couldn't identify."

    # Function to compute the carbon chain length
    def compute_carbon_chain(atom, visited=None):
        if visited is None:
            visited = set()
        max_length = 0
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in visited and neighbor.GetSymbol() == 'C':
                visited.add(neighbor.GetIdx())
                chain_length = 1 + compute_carbon_chain(neighbor, visited)
                max_length = max(max_length, chain_length)
                visited.remove(neighbor.GetIdx())
        return max_length

    # Initiate chain length check from the identified carbon
    carbon_chain_length = compute_carbon_chain(thioester_carbon)

    # Validate the length using the defined criteria
    if carbon_chain_length > 22:
        return True, f"Fatty acyl chain length is {carbon_chain_length} (> 22 carbons, valid)"
    else:
        return False, f"Fatty acyl chain length is {carbon_chain_length} (<= 22 carbons, not valid)"