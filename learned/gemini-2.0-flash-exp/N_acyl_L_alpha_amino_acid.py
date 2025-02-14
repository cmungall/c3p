"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for L-alpha-amino acid pattern.
    # Note: [C@H] denotes L stereochemistry.
    amino_acid_pattern = Chem.MolFromSmarts("[C@H](C(=O)O)N")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Molecule does not contain L-alpha-amino acid core"

    # Check for N-acyl group attached to the alpha amino
    acyl_group_pattern = Chem.MolFromSmarts("[NX3][CX3]=O")

    
    matches = mol.GetSubstructMatches(acyl_group_pattern)
    
    
    if not matches:
        return False, "Molecule does not contain an N-acyl group"

    
    # Find the specific N-H of the amino group.
    amino_matches = mol.GetSubstructMatches(amino_acid_pattern)

    if not amino_matches:
        return False, "Error: no amino acid group found even after initial match"

    # Get the atom index of the N in the amino group.
    amino_n_index = amino_matches[0][2]

    # check that the N of the amide group is the same as that of the amino group.
    amide_found = False

    for match in matches:
      for n_index in match:
        # check for atom type
        atom = mol.GetAtomWithIdx(n_index)
        if atom.GetAtomicNum() == 7 and n_index == amino_n_index:
          amide_found = True
          break
      if amide_found:
          break
    if not amide_found:
      return False, "N-acyl group is not attached to the alpha-amino group"


    return True, "Molecule is an N-acyl-L-alpha-amino acid"