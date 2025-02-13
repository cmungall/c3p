"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
"""
Classifies: CHEBI:33566 unsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    An unsaturated fatty acyl-CoA results from the formal condensation of the thiol group
    of coenzyme A with the carboxy group of any unsaturated fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA backbone pattern
    coa_pattern = Chem.MolFromSmarts("C1OC(COP(O)(=O)OP(O)(=O)OCP(O)(=O)OC2OC(N3C=NC4=C3N=CN=C4)C(O)C2OP(O)(O)=O)OC1COP(O)(O)=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA backbone found"
    
    # Look for ester bond (-C(=O)-S-) connecting fatty acid and CoA
    ester_pattern = Chem.MolFromSmarts("C(=O)SC")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester bond connecting fatty acid and CoA"
    
    # Look for carbon chain (at least 4 carbons) with at least one C=C bond
    chain_pattern = Chem.MolFromSmarts("C~C~C~C~C~C=C~C")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No unsaturated carbon chain (C=C) found"
    
    # Count rotatable bonds to verify long enough chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Carbon chain too short for fatty acid"
    
    # Count carbons and check molecular weight range
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if c_count < 12 or mol_wt < 500:
        return False, "Molecular weight or carbon count too low for fatty acyl-CoA"
    
    return True, "Contains unsaturated fatty acid chain linked to CoA backbone via ester bond"