"""
Classifies: CHEBI:37141 organobromine compound
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str) -> tuple[bool | None, str | None]:
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound contains at least one carbon-bromine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple[bool, str]: A tuple containing a boolean indicating if it is an organobromine compound
                          and a string describing the reason.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
         return None, None

    # Check for the presence of a C-Br bond using SMARTS.
    # Match any Bromine atom that's not attached to Hydrogen or another halogen
    carbon_bromine_pattern = Chem.MolFromSmarts("[Br][!H,F,Cl,I,At]")
    if mol.HasSubstructMatch(carbon_bromine_pattern):
        return True, "Contains at least one carbon-bromine bond"
    else:
        return False, "Does not contain a carbon-bromine bond"