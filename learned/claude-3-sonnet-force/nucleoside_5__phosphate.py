"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: CHEBI:37563 nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate is a ribosyl or deoxyribosyl derivative of a pyrimidine or purine base
    in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for pyrimidine or purine base
    base_pattern = Chem.MolFromSmarts("[*]1[*]c2[nH]c[nH]c2[*]c1[*]") # Purine
    base_pattern = base_pattern.Union(Chem.MolFromSmarts("[*]1[*]c[nH]c([nH]c1[*])[*]")) # Pyrimidine
    if not mol.HasSubstructMatch(base_pattern):
        return False, "No purine or pyrimidine base found"
    
    # Look for ribose/deoxyribose ring
    ribose_pattern = Chem.MolFromSmarts("[OX2]C1C(C(C(O[CX4H2]([OX2H0])(=[OX1])O)O1)O)O") # Ribose
    deoxyribose_pattern = Chem.MolFromSmarts("[OX2]C1C(C(C(O[CX4H2]([OX2H0])(=[OX1])O)O1)O)N") # Deoxyribose
    if not (mol.HasSubstructMatch(ribose_pattern) or mol.HasSubstructMatch(deoxyribose_pattern)):
        return False, "No ribose or deoxyribose ring found"
    
    # Look for phosphate group(s) attached to C-5 of ribose ring
    phosphate_pattern = Chem.MolFromSmarts("[OX2H0]P(=[OX1])(O)(O)")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"
    
    # Check if phosphate group is attached to C-5 of ribose ring
    for match in phosphate_matches:
        ribose_atoms = mol.GetAtoms()[match[0]].GetNeighbors()
        for atom in ribose_atoms:
            if atom.GetAtomicNum() == 8 and atom.HasBondWithOrder(6, 1):
                ribose_atoms = atom.GetNeighbors()
                for ribose_atom in ribose_atoms:
                    if ribose_atom.GetAtomicNum() == 8 and ribose_atom.HasBondWithOrder(6, 1):
                        ribose_atoms = ribose_atom.GetNeighbors()
                        for ribose_atom in ribose_atoms:
                            if ribose_atom.GetAtomicNum() == 8 and ribose_atom.HasBondWithOrder(6, 1):
                                ribose_atoms = ribose_atom.GetNeighbors()
                                for ribose_atom in ribose_atoms:
                                    if ribose_atom.GetAtomicNum() == 8 and ribose_atom.HasBondWithOrder(6, 1):
                                        ribose_atoms = ribose_atom.GetNeighbors()
                                        for ribose_atom in ribose_atoms:
                                            if ribose_atom.GetAtomicNum() == 6 and ribose_atom.GetDegree() == 4:
                                                return True, "Molecule contains a purine or pyrimidine base attached to a ribose or deoxyribose ring with a phosphate group at C-5"
    
    return False, "Phosphate group not attached to C-5 of ribose ring"