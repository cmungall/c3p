"""
Classifies: CHEBI:17297 UDP-sugar
"""
from rdkit import Chem

def is_UDP_sugar(smiles: str):
    """
    Determine if a molecule is a UDP-sugar based on its SMILES string. A UDP-sugar is a pyrimidine nucleotide-sugar with UDP as the nucleotide component attached via an anomeric diphosphate linkage.

    Args:
        smiles (str): The SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a UDP-sugar, False otherwise.
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify the uridine moiety pattern: uracil ring and ribose
    # Using a generalized uridine pattern
    uridine_pattern = Chem.MolFromSmarts("n1cc(=O)[nH]c(=O)c1C2OC(C(O)C(O)C2)CO")
    if not mol.HasSubstructMatch(uridine_pattern):
        return False, "Uridine moiety not found"

    # Identify the diphosphate linkage pattern
    diphosphate_pattern = Chem.MolFromSmarts("O=P(O)([O-])OP(O)([O-])=O")
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "Diphosphate linkage not found"

    return True, "Contains UDP moiety with a sugar via diphosphate linkage"