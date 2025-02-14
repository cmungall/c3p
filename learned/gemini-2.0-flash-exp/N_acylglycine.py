"""
Classifies: CHEBI:16180 N-acylglycine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine is a glycine where the nitrogen atom is acylated with another group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an N-acylglycine, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core glycine pattern with acylation on the N
    acyl_glycine_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=[OX1])[OX2]")
    
    if not mol.HasSubstructMatch(acyl_glycine_pattern):
            return False, "Molecule does not contain the N-acylglycine core structure"

    glycine_pattern = Chem.MolFromSmarts("NCC(=O)O")
    glycine_matches = mol.GetSubstructMatches(glycine_pattern)
    
    if len(glycine_matches) != 1:
        return False, f"Molecule must have exactly one glycine moiety, but found {len(glycine_matches)}"

    # Check that the N of the glycine is only connected to a carbonyl and not another amino acid
    n_atom_idx = mol.GetSubstructMatch(glycine_pattern)[0]
    n_atom = mol.GetAtomWithIdx(n_atom_idx)
    
    for neighbor in n_atom.GetNeighbors():
        if neighbor.GetSymbol() == 'C' and neighbor.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
           c_neighbors = list(neighbor.GetNeighbors())
           if any(x.GetSymbol() == 'O' for x in c_neighbors) and any(x.GetSymbol() == 'O' for x in c_neighbors): #checks for O=C
              continue # it's a carbonyl, which is correct.
           else:
              return False, "Nitrogen of glycine has an unexpected neighbor" # if no carbonyl detected
        elif neighbor.GetSymbol() == 'H':
            continue #that's okay
        else:
           return False, "Nitrogen of glycine has an unexpected neighbor"

    return True, "Molecule contains a glycine with the nitrogen acylated (N-acylglycine)"