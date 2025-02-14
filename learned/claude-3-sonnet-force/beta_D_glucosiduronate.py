"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: CHEBI:36937 beta-D-glucosiduronate
"""
from rdkit import Chem
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
    
    # Look for glucuronide substructure (glucose with carboxylate group)
    glucuronide_pattern = Chem.MolFromSmarts("[C@@H]1([C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)C([O-])=O")
    if not mol.HasSubstructMatch(glucuronide_pattern):
        return False, "No glucuronide substructure found"
    
    # Check for ring strain (glucuronides typically have pyranose ring)
    ring_atoms = mol.GetRingInfo().AtomRings()[0]
    ring_bond_orders = [mol.GetBondBetweenAtoms(ring_atoms[i], ring_atoms[i-1]).GetBondType() 
                        for i in range(len(ring_atoms))]
    if 1.0 not in ring_bond_orders:  # No single bonds in ring
        return False, "Abnormal ring structure for glucuronide"
    
    # Count hydroxy groups (should be 4)
    num_hydroxy = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1)
    if num_hydroxy != 4:
        return False, f"Expected 4 hydroxy groups, got {num_hydroxy}"
    
    # Check molecular weight (typically 200-400 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 400:
        return False, "Molecular weight outside typical range for glucuronide"
    
    return True, "Contains glucuronide substructure (glucose with carboxylate group)"