"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: CHEBI:82198 wax
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    A wax is defined as an organic compound or mixture of compounds
    that is composed of long-chain molecules and is malleable at ambient temperatures.
    Commonly, waxes are esters formed from long-chain fatty acids and long-chain alcohols.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    
    # Attempt to parse SMILES with error handling
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception as e:
        return False, f"SMILES parsing error: {str(e)}"
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify ester groups using SMARTS pattern
    ester_smarts = '[CX3](=O)[OX2H0][#6]'
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # Define threshold for long-chain (number of carbons)
    chain_length_threshold = 12

    # For each ester group, split the molecule and analyze fragments
    for match in ester_matches:
        carbonyl_c_idx = match[0]
        ester_o_idx = match[2]

        # Break the ester bond (between carbonyl carbon and oxygen)
        bonds_to_break = [(carbonyl_c_idx, ester_o_idx)]
        frags = Chem.FragmentOnBonds(mol, [mol.GetBondBetweenAtoms(*pair).GetIdx() for pair in bonds_to_break])
        frag_mols = Chem.GetMolFrags(frags, asMols=True, sanitizeFrags=True)

        # Initialize counts for long chains
        long_chain_count = 0

        # Analyze each fragment
        for frag in frag_mols:
            carbon_count = count_carbons(frag)
            if carbon_count >= chain_length_threshold:
                long_chain_count += 1

        # If both fragments have long aliphatic chains, it's a wax
        if long_chain_count >= 2:
            return True, f"Ester with two long carbon chains found (long chains: {long_chain_count})"

    return False, "No ester with two long carbon chains found"

def count_carbons(mol):
    """
    Counts the number of carbon atoms in a molecule.

    Args:
        mol (rdkit.Chem.Mol): Molecule to analyze

    Returns:
        int: Number of carbon atoms in the molecule
    """
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    return carbon_count