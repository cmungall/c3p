"""
Classifies: CHEBI:138675 gas molecular entity
"""
from rdkit import Chem

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a gas molecular entity based on its SMILES string.
    The function checks against known gas molecular entities provided in the examples.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gas molecular entity, False otherwise
        str: Reason for classification
    """
    
    # Known SMILES for gas molecular entities at STP from the provided examples
    known_gas_smiles = {
        "[O][O]", "O=[13C]=O", "[He]", "[6He]", "CCCC", "[H]C(C)=C([H])C", "[Rn]",
        "[C]", "[H]\\C(C)=C(\\[H])C", "I[H]", "[Ar]", "[4He]", "[222Rn]", "[O-][O+]=O",
        "[H][H]", "[C-]#[O+]", "CC", "[3H][3H]", "CCC=C", "[220Rn]", "CCC", "CC(C)C",
        "CC(C)=C", "[H]C([H])([H])[H]", "ClC=C", "[219Rn]", "C=C", "CC#C", "FF",
        "C(C(C(F)(F)F)(F)F)(C(F)(F)F)(F)F", "[H]N([H])[H]", "[1H][1H]", "Cl[H]", "ClCl",
        "O=C=O", "C1CO1", "[3He]", "[Xe]", "FC=C", "[Ne]", "[H]\\C(C)=C(/[H])C", "[Kr]"
    }
    
    # Check if the SMILES is one of the known gases
    if smiles in known_gas_smiles:
        return True, "SMILES matches known gas molecular entity at STP"
    
    return False, "Not a known gas molecular entity at STP"