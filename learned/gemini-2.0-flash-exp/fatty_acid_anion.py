"""
Classifies: CHEBI:28868 fatty acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is the deprotonated form of a fatty acid, characterized by a carboxylate group (-COO-) at the end of a long hydrocarbon chain (>= 6 carbons).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid anion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a terminal carboxylate group using SMARTS pattern
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "Molecule does not contain a terminal carboxylate group."

    # Find a carbon atom connected to the carboxylate group
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if not carboxylate_matches:
        return False, "Molecule does not contain a terminal carboxylate group."

    #Check the chain length to see if there are at least 6 carbons using a more general SMARTS pattern
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    carbon_chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(carbon_chain_matches) == 0:
        return False, "Carboxylate group not attached to a carbon chain with length >= 6 carbons"

    return True, "Contains a carboxylate group at the end of a carbon chain (>=6 carbons), consistent with a fatty acid anion."