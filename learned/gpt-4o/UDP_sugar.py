"""
Classifies: CHEBI:17297 UDP-sugar
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

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

    # Identify the uridine part (a pyrimidin-2-one with ribose)
    uridine_pattern = Chem.MolFromSmarts("n1c(=O)[nH]c2c1n(C[C@H]3O[C@H](COP(O)(=O)OP(O)(=O)o)C3O)c2=O")
    if not mol.HasSubstructMatch(uridine_pattern):
        return False, "Uridine moiety not found"

    # Identify the sugar part attached via diphosphate linkage
    diphosphate_sugar_pattern = Chem.MolFromSmarts("COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@@H]([C@H](O)[C@H]1O)n2ccc(=O)[nH]c2=O)[C@H](O)")
    if not mol.HasSubstructMatch(diphosphate_sugar_pattern):
        return False, "Diphosphate linkage to sugar not found"

    # Identify the anomeric carbon connection
    anomeric_linkage_pattern = Chem.MolFromSmarts("[n!H0]c1n(C[C@H]2O[C@H](CO)C(O)C(O)C2)C(NC2=O)=O1")
    if not mol.HasSubstructMatch(anomeric_linkage_pattern):
        return False, "Anomeric linkage not found with sugar"

    return True, "Contains UDP moiety with sugar via an anomeric diphosphate linkage"