"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    A 17alpha-hydroxy steroid has a hydroxyl group at carbon 17 with alpha stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the basic steroid core with a 17-hydroxy substituent (no stereochemistry requirement yet)
    # Steroid core, C17 connected to C and O
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C][C][C]2[C]1[C][C][C]3[C]2[C][C][C]4[C]3[C](O)[C][C]4")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core with a 17-hydroxy group found."
    
    #Find the C17 atom
    matches = mol.GetSubstructMatches(steroid_core_pattern)
    
    #If no matches for the substructure, return False
    if not matches:
        return False, "No steroid core match found"
    
    # We take the first match, as any match must have the correct C-17 atom with attached -OH, and only consider the first match
    match = matches[0]
    
    #Find the index of C17 atom in the match - it is the penultimate atom in the pattern
    c17_index_in_match = 15 # from the SMARTS string we used
    
    #Get atom index of c17 in the overall molecule by indexing the match object with the index
    c17_atom_index = match[c17_index_in_match]
    
    # Get the C17 atom object in the molecule
    c17_atom = mol.GetAtomWithIdx(c17_atom_index)
    
    # Get its stereo information
    c17_stereo = c17_atom.GetChiralTag()
    
    #Check if atom is chiral
    if c17_stereo == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
        return False, "C17 has no stereochemistry information available"
    
    #Check for alpha-configuration by seeing if it is defined as anti-clockwise
    if c17_stereo != Chem.rdchem.ChiralType.CHI_1:
        return False, "C17 stereochemistry not alpha"
    
    
    return True, "Contains a steroid core with a 17-alpha hydroxyl group"