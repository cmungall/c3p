"""
Classifies: CHEBI:17297 UDP-sugar
"""
"""
Classifies: CHEBI:27599 UDP-sugar

A pyrimidine nucleotide-sugar having UDP as the nucleotide component attached to an unspecified sugar via an anomeric diphosphate linkage.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a UDP-sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for UDP pattern
    udp_pattern = Chem.MolFromSmarts("C1OC(C(O)C1*N1C=NC2=C1N=CN2)COP(O)(=O)OP(O)(=O)O")
    if not mol.HasSubstructMatch(udp_pattern):
        return False, "No UDP substructure found"
    
    # Look for sugar pattern (ring with multiple oxygens)
    sugar_pattern = Chem.MolFromSmarts("[OX2r]~[CX4r]~[OX2r]")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar substructure found"
    
    # Look for anomeric diphosphate linkage
    diphosphate_pattern = Chem.MolFromSmarts("[OX2r][P]([OX2r])([OX2r])=O")
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)
    if not diphosphate_matches:
        return False, "No diphosphate linkage found"
    
    # Check if the diphosphate connects UDP and sugar
    for match in diphosphate_matches:
        udp_atom = mol.GetAtomWithIdx(match[2])
        sugar_atom = mol.GetAtomWithIdx(match[0])
        if udp_atom.IsInRingSize(6) and sugar_atom.IsInRing():
            return True, "Contains UDP linked to a sugar via an anomeric diphosphate linkage"
    
    return False, "UDP and sugar not connected via a diphosphate linkage"