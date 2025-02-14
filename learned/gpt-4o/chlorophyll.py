"""
Classifies: CHEBI:28966 chlorophyll
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    Chlorophylls are magnesium porphyrins defined by the presence of a fifth ring
    and often include long hydrophobic chains or phytol chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) True if molecule is a chlorophyll, False otherwise,
               including a reason for the classification
    """

    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for magnesium porphyrin core with a fifth ring
    mg_porphyrin_pattern = Chem.MolFromSmarts("[Mg]1:n2c(nc(nc2c3c1nccc3))(=*)*")  # Simple representation

    # Check for the presence of the magnesium porphyrin core
    if not mol.HasSubstructMatch(mg_porphyrin_pattern):
        return False, "No magnesium porphyrin core with an additional ring found"

    # Check for ester or related hydrophobic chain pattern (simplified pattern)
    hydrophobic_chain_pattern = Chem.MolFromSmarts("CCCCCCCC")  # Represents a long alkyl chain
    if not mol.HasSubstructMatch(hydrophobic_chain_pattern):
        return False, "No long hydrophobic chain detected"

    return True, "Contains magnesium porphyrin core with fifth ring and hydrophobic chain"