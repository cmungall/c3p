"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    A medium-chain fatty acyl-CoA(4-) is an acyl-CoA molecule with an aliphatic acyl chain of 6-12 carbons,
    containing no rings, no aromatic systems, and composed only of C, H, and O atoms in the acyl chain.
    The phosphate and diphosphate groups are deprotonated (negative charges).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphorus atoms (CoA derivative)
    contains_phosphorus = any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms())
    if not contains_phosphorus:
        return False, "No phosphorus atoms found; molecule is not a CoA derivative"

    # Check for adenine moiety (part of CoA)
    adenine_pattern = Chem.MolFromSmarts('n1cnc2c1ncnc2N')
    adenine_matches = mol.GetSubstructMatches(adenine_pattern)
    if not adenine_matches:
        return False, "No adenine moiety found; molecule is not a CoA derivative"

    # Look for thioester linkage C(=O)-S
    thioester_pattern = Chem.MolFromSmarts('C(=O)[S;D2]')  # Ensure sulfur is divalent
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Use the first thioester match
    acyl_carbon_idx = thioester_matches[0][0]  # Carbonyl carbon
    sulfur_idx = thioester_matches[0][2]       # Sulfur atom

    # Break bond between acyl carbon and sulfur to isolate acyl chain
    bond = mol.GetBondBetweenAtoms(acyl_carbon_idx, sulfur_idx)
    if bond is None:
        return False, "Cannot find bond between acyl carbon and sulfur"

    bond_idx = bond.GetIdx()
    mol_frag = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=False)

    # Get fragments
    frags = Chem.GetMolFrags(mol_frag, asMols=True, sanitizeFrags=True)
    acyl_chain_frag = None
    for frag in frags:
        # Fragment that does not contain phosphorus is the acyl chain
        has_phosphorus = any(atom.GetAtomicNum() == 15 for atom in frag.GetAtoms())
        if not has_phosphorus:
            acyl_chain_frag = frag
            break

    if acyl_chain_frag is None:
        return False, "Could not isolate acyl chain fragment"

    # Check for rings in acyl chain
    ring_info = acyl_chain_frag.GetRingInfo()
    if ring_info.NumRings() > 0:
        return False, "Acyl chain contains rings, not a fatty acyl chain"

    # Check for aromatic atoms in acyl chain
    if any(atom.GetIsAromatic() for atom in acyl_chain_frag.GetAtoms()):
        return False, "Acyl chain contains aromatic atoms, not a fatty acyl chain"

    # Check for allowed atoms (C, H, O) in acyl chain
    allowed_atoms = {6, 1, 8}  # Carbon, Hydrogen, Oxygen
    for atom in acyl_chain_frag.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Acyl chain contains disallowed atom: {atom.GetSymbol()}"

    # Count carbons in acyl chain (excluding carbonyl carbon already counted)
    carbon_count = sum(1 for atom in acyl_chain_frag.GetAtoms() if atom.GetAtomicNum() == 6)

    if 6 <= carbon_count <= 12:
        return True, f"Acyl chain has {carbon_count} carbons and is aliphatic; within medium-chain length (6-12)"
    else:
        return False, f"Acyl chain has {carbon_count} carbons, not within medium-chain length (6-12)"