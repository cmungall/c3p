"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
"""
Classifies: CHEBI:57499 2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA is an unsaturated fatty acyl-CoA in which the S-acyl group contains
    a double bond between positions 2 and 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA backbone pattern
    coa_pattern = Chem.MolFromSmarts("[C@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O))=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA backbone found"
    
    # Look for 2-enoyl pattern (C=C-C(=O))
    enoyl_pattern = Chem.MolFromSmarts("C=CC(=O)")
    enoyl_matches = mol.GetSubstructMatches(enoyl_pattern)
    if not enoyl_matches:
        return False, "No 2-enoyl group found"
    
    # Check that the double bond is between positions 2 and 3
    for match in enoyl_matches:
        atom1 = mol.GetAtomWithIdx(match[0])
        atom2 = mol.GetAtomWithIdx(match[1])
        atom3 = mol.GetAtomWithIdx(match[2])
        
        if atom1.GetSymbol() == 'C' and atom2.GetSymbol() == 'C' and atom3.GetSymbol() == 'C':
            bond1 = mol.GetBondBetweenAtoms(match[0], match[1])
            bond2 = mol.GetBondBetweenAtoms(match[1], match[2])
            
            if bond1.GetBondType() == Chem.BondType.DOUBLE and bond2.GetBondType() == Chem.BondType.SINGLE:
                return True, "Molecule contains a 2-enoyl-CoA group"
    
    return False, "Double bond not between positions 2 and 3"