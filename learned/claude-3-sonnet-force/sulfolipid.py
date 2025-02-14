"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: CHEBI:35967 sulfolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid is defined as a compound containing a sulfonic acid residue
    joined by a carbon-sulfur bond to a lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sulfonic acid group (-S(=O)(=O)O)
    sulfonic_acid_pattern = Chem.MolFromSmarts("S(=O)(=O)O")
    sulfonic_acid_matches = mol.GetSubstructMatches(sulfonic_acid_pattern)
    if not sulfonic_acid_matches:
        return False, "No sulfonic acid group found"
    
    # Look for carbon-sulfur bond (C-S)
    c_s_bond_pattern = Chem.MolFromSmarts("[C]-[S]")
    c_s_bond_matches = mol.GetSubstructMatches(c_s_bond_pattern)
    if not c_s_bond_matches:
        return False, "No carbon-sulfur bond found"
    
    # Check if sulfonic acid group is attached to lipid chain
    lipid_chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]")
    lipid_chain_matches = mol.GetSubstructMatches(lipid_chain_pattern)
    if not lipid_chain_matches:
        return False, "No lipid chain found"
    
    # Check if sulfonic acid and lipid chain are connected via carbon-sulfur bond
    for s_idx in sulfonic_acid_matches:
        for c_idx in c_s_bond_matches:
            if s_idx in c_idx:
                sulfur_atom = mol.GetAtomWithIdx(s_idx)
                for bond in sulfur_atom.GetBonds():
                    if bond.GetOtherAtom(sulfur_atom).GetAtomicNum() == 6:  # Carbon
                        carbon_atom = bond.GetOtherAtom(sulfur_atom)
                        is_lipid_connected = False
                        for lipid_match in lipid_chain_matches:
                            if carbon_atom.GetIdx() in lipid_match:
                                is_lipid_connected = True
                                break
                        if is_lipid_connected:
                            return True, "Contains sulfonic acid group joined by a carbon-sulfur bond to a lipid chain"
    
    return False, "Sulfonic acid group not connected to lipid chain via carbon-sulfur bond"