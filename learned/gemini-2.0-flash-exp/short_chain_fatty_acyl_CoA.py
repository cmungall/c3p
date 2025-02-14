"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    A short-chain fatty acyl-CoA is a CoA molecule linked to a short-chain fatty acid via a thioester.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core CoA substructure (pantetheine part)
    coa_core_smarts = "[SX2][CH2X4][CH2X4][NX3][CX3](=[OX1])[NX3][CH2X4][CH2X4][CX3](=[OX1])[C@H]([OX2])[CX4]([CH3X4])([CH3X4])[CX4][OX2][PX4]"
    coa_core_pattern = Chem.MolFromSmarts(coa_core_smarts)
    if not mol.HasSubstructMatch(coa_core_pattern):
         return False, "CoA core not found"

    # Define the thioester bond
    thioester_smarts = "[CX3](=[OX1])[SX2]"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester bond found"
    
    # Check if thioester is attached to CoA
    coa_match = mol.GetSubstructMatch(coa_core_pattern)
    attached_to_coa = False
    for thioester_match in thioester_matches:
        for atom_index in thioester_match:
            atom = mol.GetAtomWithIdx(atom_index)
            if atom.GetSymbol() == "S":
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() in coa_match:
                       attached_to_coa = True
                       break
                if attached_to_coa:
                    break
        if attached_to_coa:
            break
    
    if not attached_to_coa:
        return False, "Thioester not attached to CoA core"

    # Find the acyl chain and count carbons.  We do this by looking at neighbors of the carbonyl carbon
    # Note that this does not guarantee that it is a single chain, just a chain
    acyl_carbon_count = 0
    for thioester_match in thioester_matches:
        carbonyl_carbon_index = -1
        for atom_index in thioester_match:
            atom = mol.GetAtomWithIdx(atom_index)
            if atom.GetSymbol() == "C" and atom.GetHybridization() == Chem.HybridizationType.SP2:
                carbonyl_carbon_index = atom_index
                break

        if carbonyl_carbon_index != -1:
             carbonyl_atom = mol.GetAtomWithIdx(carbonyl_carbon_index)
             chain_atoms = set()
             queue = [carbonyl_carbon_index]

             while queue:
                 current_index = queue.pop(0)
                 if current_index in chain_atoms:
                     continue
                 chain_atoms.add(current_index)
                 current_atom = mol.GetAtomWithIdx(current_index)
                 if current_atom.GetSymbol() == "C" :
                    acyl_carbon_count += 1
                    for neighbor in current_atom.GetNeighbors():
                        if neighbor.GetSymbol() == "C" and neighbor.GetIdx() not in chain_atoms and neighbor.GetIdx() != carbonyl_carbon_index : # avoid going back to carbonyl carbon
                            queue.append(neighbor.GetIdx())
                
    # Remove the carbonyl carbon itself
    acyl_carbon_count -= 1

    if acyl_carbon_count < 1 or acyl_carbon_count > 7:
        return False, f"Fatty acid chain has {acyl_carbon_count} carbons, must be between 1 and 7"
    
    return True, "Molecule is a short-chain fatty acyl-CoA"