"""
Classifies: CHEBI:4194 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a monosaccharide with six carbons and the D configuration at carbon 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count != 6:
        return False, f"Not a hexose. Found {c_count} carbon atoms"
    if o_count < 5:
        return False, f"Not a hexose. Found {o_count} oxygen atoms"


    # Identify carbons and chiral center
    
    # Look for a primary alcohol as a common endgroup
    primary_alcohol_pattern = Chem.MolFromSmarts("[CH2X4]-O")
    alcohol_matches = mol.GetSubstructMatches(primary_alcohol_pattern)
    
    if len(alcohol_matches) != 1:
         return False, "Not a hexose, no primary alcohol found"
    
    
    # Locate the primary alcohol carbon and label it as C6.
    c6_index = alcohol_matches[0][0]
    
    # Find the carbon atoms connected to C6
    neighbors_c6 = [atom.GetIdx() for atom in mol.GetAtomWithIdx(c6_index).GetNeighbors()]
    c5_index = -1
    for neighbor in neighbors_c6:
        if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6 and neighbor != c6_index:
             c5_index = neighbor
             break
    
    if c5_index == -1:
        return False, "Cannot identify carbon C5."

    # Find the other four carbons by looking for adjacent carbons on the chain
    c4_index = -1
    c3_index = -1
    c2_index = -1
    c1_index = -1
    
    neighbors_c5 = [atom.GetIdx() for atom in mol.GetAtomWithIdx(c5_index).GetNeighbors()]
    for neighbor in neighbors_c5:
        if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6 and neighbor != c6_index:
            c4_index = neighbor
            break
    if c4_index == -1:
         return False, "Cannot identify carbon C4."
        
    neighbors_c4 = [atom.GetIdx() for atom in mol.GetAtomWithIdx(c4_index).GetNeighbors()]
    for neighbor in neighbors_c4:
        if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6 and neighbor != c5_index:
            c3_index = neighbor
            break
    if c3_index == -1:
         return False, "Cannot identify carbon C3."
    
    neighbors_c3 = [atom.GetIdx() for atom in mol.GetAtomWithIdx(c3_index).GetNeighbors()]
    for neighbor in neighbors_c3:
        if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6 and neighbor != c4_index:
            c2_index = neighbor
            break
    if c2_index == -1:
         return False, "Cannot identify carbon C2."
        
    neighbors_c2 = [atom.GetIdx() for atom in mol.GetAtomWithIdx(c2_index).GetNeighbors()]
    for neighbor in neighbors_c2:
        if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6 and neighbor != c3_index:
             c1_index = neighbor
             break
    if c1_index == -1:
        return False, "Cannot identify carbon C1."

    # Get chirality at C5
    c5_chirality = mol.GetAtomWithIdx(c5_index).GetChiralTag()

    # Check if chiral center has the R configuration, use absolute chirality to handle cases where chirality can be S or R
    if c5_chirality == Chem.ChiralType.CHI_TETRAHEDRAL_CW or c5_chirality == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
         # R configuration of the chiral center is encoded by absolute chirality
         # In RDKit, "clockwise" is R when looking along the bond from chiral center to substituent with lowest priority.
         # In carbohydrates, this is the substituent with the highest atomic number.
         
         # Obtain all four neighbors
         neighbors_c5 = [atom.GetIdx() for atom in mol.GetAtomWithIdx(c5_index).GetNeighbors()]
         
         # Create a list of neighbors of c5 sorted by atomic number (highest first)
         neighbor_atoms = sorted([(neighbor, mol.GetAtomWithIdx(neighbor).GetAtomicNum()) for neighbor in neighbors_c5], key=lambda x: x[1], reverse=True)

         # Find the atom with the lowest priority and place it at the end of the list.
         lowest_priority = [neighbor[0] for neighbor in neighbor_atoms if neighbor[1] == 1]
         if len(lowest_priority) == 0:
             lowest_priority = [neighbor_atoms[-1][0]]

         if len(lowest_priority) == 1:
             # Find the three other neighbors and add them to the list
             neighbor_list = [neighbor[0] for neighbor in neighbor_atoms if neighbor[0] != lowest_priority[0]]
             neighbor_list.append(lowest_priority[0])

             # Determine R or S configuration for comparison by setting the chiral tag to CHI_TETRAHEDRAL_CW
             new_mol = Chem.Mol(mol)
             new_mol.GetAtomWithIdx(c5_index).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
             
             # Compare the configurations
             if Chem.MolToSmiles(new_mol, isomericSmiles = True) == Chem.MolToSmiles(mol, isomericSmiles=True):
                   if c5_chirality != Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                     return False, "Configuration of C5 is not D"
             else:
                 if c5_chirality == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                      return False, "Configuration of C5 is not D"

         else:
              return False, "C5 has multiple hydrogens and is not chiral."
    else:
         return False, "C5 is not a chiral center."

    return True, "D-hexose configuration confirmed"