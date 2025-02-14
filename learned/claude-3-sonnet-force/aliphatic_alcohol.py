"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: CHEBI:26003 aliphatic alcohol
An alcohol derived from an aliphatic compound.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol is a molecule containing at least one aliphatic (non-aromatic)
    carbon chain and at least one hydroxyl (-OH) group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for hydroxyl group (-OH)
    if not any(atom.GetSymbol() == 'O' and 
               atom.GetHybridization() == Chem.HybridizationType.SP3 and 
               sum(mol.GetAtomWithIdx(i).GetTotalNumHs() for i in atom.GetNeighbors()) == 1
               for atom in mol.GetAtoms()):
        return False, "No hydroxyl (-OH) group found"

    # Check for aliphatic (non-aromatic) carbon chain
    if not any(bond.GetIsConjugated() is False and 
               bond.GetBeginAtom().GetSymbol() == 'C' and 
               bond.GetEndAtom().GetSymbol() == 'C'
               for bond in mol.GetBonds()):
        return False, "No aliphatic carbon chain found"

    # Passed all checks
    return True, "Contains at least one aliphatic carbon chain and one hydroxyl (-OH) group"