"""
Classifies: CHEBI:28868 fatty acid anion
"""
"""
Classifies: CHEBI:33416 fatty acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is the conjugate base of a fatty acid, arising from deprotonation
    of the carboxylic acid group of the corresponding fatty acid.

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

    # Look for carboxylate group (-COO-)
    carboxylate_pattern = Chem.MolFromSmarts("[C]([O-])=O")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"

    # Look for long carbon chain attached to the carboxylate
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if len(carbon_chain_matches) < 1:
        return False, "Missing carbon chain attached to carboxylate"

    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Carbon chain too short to be a fatty acid"

    # Check molecular weight - fatty acids typically >100 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for fatty acid anion"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 8:
        return False, "Too few carbons for fatty acid anion"
    if o_count < 2:
        return False, "Must have at least 2 oxygens (carboxylate group and other functionality)"

    return True, "Contains carboxylate group (-COO-) attached to a long carbon chain"