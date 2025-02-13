"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
"""
Classifies: lipid hydroperoxide
Definition: Any lipid carrying one or more hydroperoxy substituents.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_match, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for hydroperoxy group (-OO)
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2][OX2H]")
    hydroperoxy_matches = mol.GetSubstructMatches(hydroperoxy_pattern)
    if not hydroperoxy_matches:
        return False, "No hydroperoxy (-OO) group found"

    # Check for lipid characteristics
    
    # Look for carboxylic acid or ester group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH,OR]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid or ester group found"
    
    # Count carbons (lipids typically have >8 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 8:
        return False, "Carbon chain too short for a lipid"
    
    # Look for carbon chain
    alkyl_chain = Chem.MolFromSmarts("[CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(alkyl_chain):
        return False, "No alkyl chain found"

    # Count hydroperoxy groups
    num_hydroperoxy = len(hydroperoxy_matches)
    
    # Check for double bonds (most lipid hydroperoxides are unsaturated)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    has_double_bonds = len(double_bond_matches) > 0
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:  # Arbitrary minimum weight for a lipid hydroperoxide
        return False, "Molecular weight too low for a lipid hydroperoxide"

    # Construct reason string
    reason_parts = []
    reason_parts.append(f"Contains {num_hydroperoxy} hydroperoxy group(s)")
    if has_double_bonds:
        reason_parts.append(f"has {len(double_bond_matches)} double bond(s)")
    reason_parts.append(f"contains {carbon_count} carbons")
    reason = ", ".join(reason_parts)
    
    return True, reason