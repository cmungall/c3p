"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies: CHEBI:35271 carbamate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is any ester of carbamic acid or its N-substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbamate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carbamate ester functional group (-O-C(=O)-N-)
    carbamate_pattern = Chem.MolFromSmarts("[OX2]C(=O)N")
    if not mol.HasSubstructMatch(carbamate_pattern):
        return False, "No carbamate ester functional group found"

    # Check for aromatic rings
    rings = mol.GetRingInfo().AtomRings()
    has_aromatic_ring = any(mol.GetAtomWithIdx(idx).GetIsAromatic() for ring in rings for idx in ring)

    # Check for alkyl groups
    has_alkyl_groups = any(atom.GetHybridization() == Chem.HybridizationType.SP3 and atom.GetSymbol() == 'C' for atom in mol.GetAtoms())

    # Carbamate esters often (but not always) have aromatic rings and/or alkyl groups
    if has_aromatic_ring or has_alkyl_groups:
        return True, "Contains carbamate ester functional group and common structural features"
    else:
        return True, "Contains carbamate ester functional group but lacks typical structural features"