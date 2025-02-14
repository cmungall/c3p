"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate is defined as at least two biological molecules covalently linked.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for common biological building blocks
    amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0][CX4][CX3](=[OX1])[OX2;H1,H0]") # Covers both N- and C- termini
    sugar_pattern = Chem.MolFromSmarts("OC[C,c]1([O,N,S][C,c]([O,N,S])[C,c]([O,N,S])[C,c]([O,N,S])[C,c]1[O,N,S])") # Covers many sugar types
    nucleotide_base_pattern = Chem.MolFromSmarts("c1[nc][nc][nc]1") # covers purine and pyrimidine bases
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([O])([O])[O]") # Covers phosphates
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") # A crude definition, but OK for now
    
    patterns = [amino_acid_pattern, sugar_pattern, nucleotide_base_pattern, phosphate_pattern, fatty_acid_pattern]

    # Find matches for each pattern
    matches = []
    for pattern in patterns:
      if mol.HasSubstructMatch(pattern):
        matches.append(pattern)

    # Check if we have at least 2 different classes of substructures.
    if len(matches) < 2:
      return False, "Does not contain at least two different bio-molecule substructures"


    # Check for covalent linkages
    # Find all substructure atom matches for each substructure
    all_substructure_atoms = []
    for pattern in matches:
      match_indices = mol.GetSubstructMatches(pattern)
      if len(match_indices) > 0: # Only process if a match is found
          for match in match_indices:
            all_substructure_atoms.append(set(match))
    
    
    # Check for atom overlap across all substructures
    # If two substructures share atoms, they are covalently linked
    linked = False
    if len(all_substructure_atoms) > 1:
        for i in range(len(all_substructure_atoms)):
            for j in range(i+1, len(all_substructure_atoms)):
                if len(all_substructure_atoms[i].intersection(all_substructure_atoms[j])) > 0:
                    linked = True
                    break
            if linked:
                break
    
    if not linked:
        return False, "Contains at least two different bio-molecule substructures, but not covalently linked."

    return True, "Contains at least two different covalently linked bio-molecule substructures."