"""
Classifies: CHEBI:20706 6-aminopurines
"""
"""
Classifies: CHEBI:26947 6-aminopurines
A 6-aminopurine is any compound having 6-aminopurine (adenine) as part of its structure.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule contains a 6-aminopurine (adenine) moiety based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a 6-aminopurine moiety, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of 6-aminopurine ring system
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    has_adenine = mol.HasSubstructMatch(adenine_pattern)
    
    if not has_adenine:
        return False, "Does not contain 6-aminopurine (adenine) moiety"
    
    # Count nitrogens to verify purine
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 5:
        return False, "Too few nitrogens for purine system"
    
    # Check for CoA or phosphate attachments (common in metabolites)
    has_coa = mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[C@@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"))
    has_phosphate = any(atom.GetDegree() == 4 and atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())
    
    # Verify mass range (optional)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 1500:
        return False, "Molecular weight out of typical range for 6-aminopurine metabolite"
    
    if has_coa or has_phosphate:
        return True, "Contains 6-aminopurine (adenine) moiety with CoA or phosphate attachment"
    else:
        return True, "Contains 6-aminopurine (adenine) moiety"