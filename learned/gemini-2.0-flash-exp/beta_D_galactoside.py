"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    A beta-D-galactoside contains a D-galactose ring with beta-configuration at its anomeric center.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-galactoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the D-galactose pyranose ring (without specifying anomeric configuration)
    galactose_ring_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@@H](O)[C@H](CO)O[C@@H]1*")
    matches = mol.GetSubstructMatches(galactose_ring_pattern)
    if not matches:
        return False, "No D-galactose ring found"
    
    # Iterate over the matches and check for beta anomeric configuration
    for match in matches:
      # The first atom in the pattern is the ring oxygen.
      ring_oxygen = match[0]
      # The carbon attached to the ring oxygen is at index 1 in the pattern
      anomeric_carbon = match[1]
      
      # Get the RDKit atom object from the index
      anomeric_carbon_atom = mol.GetAtomWithIdx(anomeric_carbon)

      #Get the chiral tag for the anomeric carbon atom. If the carbon is not chiral it returns a non-chiral tag.
      chiral_tag = anomeric_carbon_atom.GetChiralTag()

      if chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
          # It is beta
          return True, "Contains a beta-D-galactoside moiety"
          
    return False, "Not a beta-D-galactoside, based on anomeric carbon stereochemistry"