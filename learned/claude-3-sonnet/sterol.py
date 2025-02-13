"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: CHEBI:18033 sterol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is defined as any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol
    (additional carbon atoms may be present in the side chain).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 3-hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[OH1]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) != 1:
        return False, "Must have exactly one hydroxy group"
    if not mol.GetAtomWithIdx(hydroxy_matches[0][0]).IsInRing():
        return False, "Hydroxy group must be in a ring system"

    # Check for steroid backbone pattern
    steroid_pattern = Chem.MolFromSmarts("[C@]12[C@H]3[C@@H](CCC4[C@@]1(CCC5[C@@]2(CCC(=O)C5)C)C)C6=CC=CC6=C3")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Does not contain steroid backbone"

    # Check for additional side chains
    side_chain_pattern = Chem.MolFromSmarts("[C@H](CC[C@@H](C)C)C")
    side_chain_matches = mol.GetSubstructMatches(side_chain_pattern)
    if len(side_chain_matches) > 1:
        return False, "More than one side chain present"

    # Check molecular weight - sterols typically 300-500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 500:
        return False, "Molecular weight outside typical range for sterols"

    # Count rings and carbons
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_rings < 3 or c_count < 20:
        return False, "Insufficient rings or carbons for sterol"

    return True, "Contains steroid backbone with 3-hydroxy group and permitted side chains"