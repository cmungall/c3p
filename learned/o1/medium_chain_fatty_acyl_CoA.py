"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:60903 medium-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.
    A medium-chain fatty acyl-CoA is a fatty acyl-CoA resulting from the condensation of 
    coenzyme A with a medium-chain fatty acid (6 to 12 carbons).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define thioester pattern: carbonyl carbon attached to sulfur
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Assume the first match is the thioester bond we're interested in
    thioester_atom_idx = thioester_matches[0][0]  # Index of carbonyl carbon

    # Find the bond between the carbonyl carbon and sulfur
    carbonyl_atom = mol.GetAtomWithIdx(thioester_atom_idx)
    sulfur_atom = None
    for neighbor in carbonyl_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 16:  # Atomic number of sulfur
            sulfur_atom = neighbor
            break
    if sulfur_atom is None:
        return False, "Thioester sulfur atom not found"

    # Get the bond between carbonyl carbon and sulfur
    thioester_bond = mol.GetBondBetweenAtoms(carbonyl_atom.GetIdx(), sulfur_atom.GetIdx())
    if thioester_bond is None:
        return False, "Thioester bond not found"

    # Fragment the molecule at the thioester bond to isolate the fatty acyl chain
    fragmented_mol = Chem.FragmentOnBonds(mol, [thioester_bond.GetIdx()])
    fragments = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)

    # Identify the fatty acyl chain fragment
    acyl_chain = None
    for frag in fragments:
        # If fragment does not contain sulfur, it is likely the acyl chain
        if not any(atom.GetAtomicNum() == 16 for atom in frag.GetAtoms()):
            acyl_chain = frag
            break
    if acyl_chain is None:
        return False, "Fatty acyl chain not found"

    # Count the number of carbons in the acyl chain
    carbon_count = sum(1 for atom in acyl_chain.GetAtoms() if atom.GetAtomicNum() == 6)

    # Medium-chain fatty acids have 6 to 12 carbons
    if 6 <= carbon_count <= 12:
        return True, f"Contains fatty acyl chain with {carbon_count} carbons (medium-chain)"
    else:
        return False, f"Fatty acyl chain with {carbon_count} carbons is not medium-chain"