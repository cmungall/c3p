"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:36345 long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import BondType

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A long-chain fatty acyl-CoA results from the condensation of CoA with a long-chain (C13-C22) fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for CoA substructure
    coa_pattern = Chem.MolFromSmarts("C(C)(COP(=O)([O-])OP(=O)([O-])OC1OC(n2cnc3c(N)ncnc23)C(O)C1OP(=O)([O-]))")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No coenzyme A substructure found"
    
    # Check for carbonyl group
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if len(carbonyl_matches) != 1:
        return False, f"Expected 1 carbonyl group, found {len(carbonyl_matches)}"
    
    # Check for thioester bond between carbonyl and CoA
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester bond between carbonyl and CoA found"
    
    # Check for long-chain fatty acid
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) == 0:
        return False, "No long fatty acid chain found"
    
    # Count carbon atoms in fatty acid chain
    fatty_acid_atoms = [mol.GetAtomWithIdx(match[0]).GetAtomicNum() for match in fatty_acid_matches]
    n_carbons = sum(1 for atom in fatty_acid_atoms if atom == 6)
    if n_carbons < 13 or n_carbons > 22:
        return False, f"Fatty acid chain has {n_carbons} carbons, expected 13-22"
    
    # Check for unsaturations/double bonds in fatty acid chain
    fatty_acid_bonds = [mol.GetBondBetweenAtoms(match[i], match[i+1]).GetBondType() for match in fatty_acid_matches for i in range(len(match)-1)]
    n_double_bonds = sum(1 for bond in fatty_acid_bonds if bond == BondType.DOUBLE)
    
    return True, f"Long-chain fatty acid with {n_carbons} carbons and {n_double_bonds} double bonds attached to CoA via thioester"