"""
Classifies: CHEBI:17984 acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is a thioester formed between coenzyme A and a carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the CoA moiety (focusing on key parts - pantothenic acid and the mercaptoethylamine)
    # This pattern might need refinement based on more examples, but seems to be ok
    # The sulfur atom is not specified in the CoA SMARTS as it is used to connect to the rest of the molecule
    coa_pattern = Chem.MolFromSmarts("CCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([OH])OP(=O)([OH])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)([OH])O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A substructure not found"

    # Define a SMARTS pattern for the thioester group (C(=O)-S)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
         return False, "Thioester group not found"

    # Check if the sulfur of the thioester is directly connected to the CoA
    
    # Iterate over the thioester matches.
    # For each match, get the sulfur atom, then iterate
    # over the atoms connected to it. If one of them matches
    # an atom from the CoA part of the molecule, we have found an acyl-CoA
    
    
    coa_match = mol.GetSubstructMatch(coa_pattern)
    is_connected = False
    for match in thioester_matches:
        thioester_sulfur_atom = mol.GetAtomWithIdx(match[1])
        for neighbor in thioester_sulfur_atom.GetNeighbors():
            if neighbor.GetIdx() in coa_match:
                is_connected = True
                break
        if is_connected:
            break
    
    if not is_connected:
        return False, "Thioester not connected to the CoA moiety"


    return True, "Molecule contains Coenzyme A and a thioester connected to a carboxylic acid."