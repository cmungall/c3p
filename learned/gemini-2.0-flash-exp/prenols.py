"""
Classifies: CHEBI:26244 prenols
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    Prenols are characterized by a repeating isoprene unit (C5H8) with a terminal alcohol group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Initial filter: must have an alcohol and at least one double bond
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No alcohol group found"
    
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    if not mol.HasSubstructMatch(double_bond_pattern):
       return False, "No C=C double bond found"
    

    # Core isoprene unit substructure:
    isoprene_core = Chem.MolFromSmarts("C[CX3]([CX4])=[CX3][CX4]") #methyl-C=C-C
    if not mol.HasSubstructMatch(isoprene_core):
       return False, "No core isoprene unit found."
    

    # Count carbons, hydrogens, and oxygens:
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    hydrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Check for sufficient carbons for the isoprene unit(s).
    if carbon_count < 5 :
      return False, "Too few carbons for an isoprene unit"

    # Check chain length based on rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
        return False, "Chain too short, less than 2 rotatable bonds"


    # Termination: verify we have only ONE alcohol group
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    if len(alcohol_matches) != 1:
         return False, f"Must have one alcohol group, found {len(alcohol_matches)}"

    # Check number of isoprene units (approximate)
    # Number of isoprene units is approximated from number of carbons. There should be at least 1, and each isoprene unit has 5 carbons.
    isoprene_units = (carbon_count - 2) / 5 #remove 2 for alcohol
    if isoprene_units < 0.8:
        return False, f"Too few isoprene units ({isoprene_units})."

    return True, "Matches prenol criteria"