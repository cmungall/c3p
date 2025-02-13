"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: CHEBI:35615 beta-D-glucosiduronate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.
    A beta-D-glucosiduronate is a carbohydrate acid derivative anion obtained by deprotonation
    of the carboxy group of any beta-D-glucosiduronic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glucuronic acid moiety ([C@H]1[C@H]([C@@H]([C@H]([C@H](O1)C([O-])=O)O)O)O)
    glucuronic_pattern = Chem.MolFromSmarts("[C@H]1[C@H]([C@@H]([C@H]([C@H](O1)C([O-])=O)O)O)O")
    if not mol.HasSubstructMatch(glucuronic_pattern):
        return False, "No glucuronic acid moiety found"
    
    # Check for beta configuration (look for cis arrangement of H and OH on ring atoms)
    beta_pattern = Chem.MolFromSmarts("[C@H]1[C@@H]([C@@H]([C@@H]([C@@H](O1)C([O-])=O)O)O)O")
    if not mol.HasSubstructMatch(beta_pattern):
        return False, "Not in beta configuration"
    
    # Check molecular weight - glucosiduronates typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for glucosiduronate"
    
    # Count carbons, oxygens, and deprotonated carboxyl groups
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    deprotonated_cooh_count = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() == -1)
    
    if c_count < 6:
        return False, "Too few carbons for glucosiduronate"
    if o_count < 6:
        return False, "Too few oxygens for glucosiduronate"
    if deprotonated_cooh_count != 1:
        return False, "Does not contain exactly one deprotonated carboxyl group"
    
    return True, "Contains glucuronic acid moiety in beta configuration with deprotonated carboxyl group"