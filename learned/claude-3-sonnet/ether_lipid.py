"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: CHEBI:38112 ether lipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    Ether lipids are lipids where one or more carbon atoms of the glycerol backbone
    are bonded to an alkyl chain via an ether linkage instead of an ester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone pattern [C-C-C]
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for ether linkage (-O-C) bonded to glycerol backbone
    ether_pattern = Chem.MolFromSmarts("[CH2X4][OX2][CX4]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if not ether_matches:
        return False, "No ether linkage found"
    
    # Check for alkyl chains (long carbon chains attached to ether linkage)
    alkyl_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    alkyl_chain_matches = []
    for match in ether_matches:
        ether_atom = match[1]  # Index of ether oxygen
        chains = mol.GetAtomWithIdx(ether_atom).GetNeighbors()
        for chain in chains:
            if chain.GetAtomicNum() == 6:  # Carbon
                chain_matches = mol.GetSubstructMatches(alkyl_chain_pattern, chain.GetIdx())
                alkyl_chain_matches.extend(chain_matches)
    
    if not alkyl_chain_matches:
        return False, "No alkyl chains found"
    
    # Check molecular weight - ether lipids typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for ether lipid"
    
    return True, "Contains glycerol backbone with one or more ether-linked alkyl chains"