"""
Classifies: CHEBI:51963 hopanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is a triterpenoid based on a hopane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the core hopane skeleton SMARTS pattern. This is 5 fused rings, with C and CC bonds.
    hopane_core_pattern = Chem.MolFromSmarts("[CX4]1[CX4]2[CX4]3[CX4]4[CX4]5[CX4]1[CX4]1[CX4]2[CX4]3[CX4]4[CX4]5")
    if not mol.HasSubstructMatch(hopane_core_pattern):
          return False, "Molecule does not contain the core hopane pentacyclic ring system."

    # Check for the presence of at least four methyl groups attached to the core ring system.
    methyl_group_pattern = Chem.MolFromSmarts("[CX4]([H])([H])([H])")
    methyl_matches = mol.GetSubstructMatches(methyl_group_pattern)
    
    if len(methyl_matches) < 4:
         return False, f"Too few methyl groups, should have at least 4, got {len(methyl_matches)}"
        
    # Check for the number of carbons (hopane is C30, but it can be higher, so at least 25 is OK)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 25:
        return False, f"Too few carbon atoms. Must have at least 25, got {num_carbons}"

    # Molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350:
      return False, "Molecular weight is too low for a hopanoid."
      
    
    return True, "Molecule contains the hopane core structure with at least 4 methyl groups and appropriate size for hopanoid"