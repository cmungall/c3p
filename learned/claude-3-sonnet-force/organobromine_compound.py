"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: CHEBI:38241 organobromine compound
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound is a compound containing at least one carbon-bromine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organobromine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule contains at least one C-Br bond
    has_c_br = False
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 35) or (a1.GetAtomicNum() == 35 and a2.GetAtomicNum() == 6):
            has_c_br = True
            break

    if has_c_br:
        return True, "Contains at least one carbon-bromine bond"
    else:
        return False, "Does not contain any carbon-bromine bonds"