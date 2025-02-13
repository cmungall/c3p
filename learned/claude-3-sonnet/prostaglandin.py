"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: CHEBI:35473 prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are naturally occurring compounds derived from the parent C20 acid, prostanoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for prostanoic acid backbone (cyclopentane ring fused to a cyclopentene ring)
    prostanoic_pattern = Chem.MolFromSmarts("[C@H]1C[C@H]([C@@H]2[C@@]1([H])CC=C2)C")
    if not mol.HasSubstructMatch(prostanoic_pattern):
        return False, "No prostanoic acid backbone found"
    
    # Look for substituents on the prostanoic acid backbone
    # Count number of substituents on the cyclopentane ring
    n_cyclopentane_subs = sum(1 for atom in mol.GetAtomWithIdx(prostanoic_pattern.GetAtomsMatchingSmarts()[0]).GetNeighbors() 
                               if atom.GetSymbol() != 'C' and atom.GetSymbol() != 'H')
    
    # Count number of substituents on the cyclopentene ring
    n_cyclopentene_subs = sum(1 for atom in mol.GetAtomWithIdx(prostanoic_pattern.GetAtomsMatchingSmarts()[1]).GetNeighbors() 
                               if atom.GetSymbol() != 'C' and atom.GetSymbol() != 'H')
    
    # Prostaglandins typically have 1-2 substituents on the cyclopentane ring and a lipid chain on the cyclopentene ring
    if n_cyclopentane_subs < 1 or n_cyclopentane_subs > 2 or n_cyclopentene_subs != 1:
        return False, "Incorrect number of substituents on prostanoic acid backbone"
    
    # Check for lipid chain (long carbon chain)
    lipid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    lipid_matches = mol.GetSubstructMatches(lipid_pattern)
    if not lipid_matches:
        return False, "No lipid chain found"
    
    # Count rotatable bonds to verify long lipid chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Lipid chain too short"
    
    # Check for common prostaglandin functional groups (hydroxy, keto, carboxyl)
    has_hydroxy = any(atom.GetSymbol() == 'O' and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 1 for atom in mol.GetAtoms())
    has_keto = any(atom.GetSymbol() == 'O' and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 2 for atom in mol.GetAtoms())
    has_carboxyl = any(atom.GetSymbol() == 'O' and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 2 and 
                       any(nbr_atom.GetSymbol() == 'O' and sum(nbr_bond.GetBondTypeAsDouble() for nbr_bond in nbr_atom.GetBonds()) == 1
                           for nbr_atom, nbr_bond in zip(atom.GetNeighbors(), atom.GetBonds()))
                       for atom in mol.GetAtoms())
    
    if not any([has_hydroxy, has_keto, has_carboxyl]):
        return False, "No common prostaglandin functional group found"
    
    return True, "Contains prostanoic acid backbone with a lipid chain and common prostaglandin functional groups"