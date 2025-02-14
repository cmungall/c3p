"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid (MUFA) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a MUFA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O,OH]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
      return False, "No carboxylic acid group found"
    
    # Count double and triple bonds
    num_double_bonds = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))
    num_triple_bonds = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C#C")))

    total_unsaturations = num_double_bonds + num_triple_bonds
    
    if total_unsaturations != 1:
        return False, f"Molecule has {total_unsaturations} double/triple bonds, should have exactly 1"
      
    # Get the number of carbons in the fatty acid chain (number of carbons - those in the carboxyl group)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    #check for rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 1:
         return False, "Contains more than one ring."
    
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if not len(ring_atoms) == 3:
            return False, "Ring contains more than 3 atoms"
            
    # Check for fatty acid chain length
    if c_count < 10 or c_count > 30 :
        return False, f"Carbon chain length ({c_count}) is not within range [10-30] for a fatty acid"
    
    # Check that every other bond is a single bond, this is complicated, but might work.
    saturated_pattern = Chem.MolFromSmarts("C-C")
    saturated_matches = len(mol.GetSubstructMatches(saturated_pattern))
    
    if num_double_bonds == 1:
      carbon_count_single_bonds = c_count -1
    elif num_triple_bonds ==1:
      carbon_count_single_bonds = c_count -2
    else:
      return False, f"No unsaturation identified"
    
    if saturated_matches != carbon_count_single_bonds -1 :
        return False, f"Chain has more than 1 unsaturated bond: {saturated_matches} bonds found, {carbon_count_single_bonds -1} expected"
            
    return True, "Monounsaturated fatty acid identified"