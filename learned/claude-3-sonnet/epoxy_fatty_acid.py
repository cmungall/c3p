"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: epoxy fatty acid
A heterocyclic fatty acid containing an epoxide ring as part of its structure.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_epoxy_fatty_acid, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[OH]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"
    
    # Check for epoxide ring
    # Three-membered ring containing oxygen
    epoxide = Chem.MolFromSmarts('[C]1[O][C]1')
    if not mol.HasSubstructMatch(epoxide):
        return False, "No epoxide ring found"
    
    # Count carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:
        return False, f"Carbon chain too short ({carbon_count} carbons) for fatty acid"
    if carbon_count > 30:
        return False, f"Carbon chain too long ({carbon_count} carbons) for fatty acid"
    
    # Count double bonds
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    # Subtract one for the carboxylic acid C=O
    alkene_count = double_bond_count - 1
    if alkene_count > 6:
        return False, f"Too many double bonds ({alkene_count}) for typical epoxy fatty acid"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 500:
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical range for epoxy fatty acids"
    
    # Count epoxide rings
    epoxide_matches = len(mol.GetSubstructMatches(epoxide))
    if epoxide_matches > 3:
        return False, f"Too many epoxide rings ({epoxide_matches})"
        
    # Check oxygen count and distribution
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 3:  # Minimum: 1 for epoxide + 2 for carboxylic acid
        return False, f"Too few oxygen atoms ({oxygen_count})"
    if oxygen_count > 10:  # Maximum: typically not more than 10
        return False, f"Too many oxygen atoms ({oxygen_count})"
    
    # Verify the presence of an aliphatic chain
    aliphatic_chain = Chem.MolFromSmarts('[CH2][CH2][CH2]')
    if not mol.HasSubstructMatch(aliphatic_chain):
        return False, "No aliphatic chain found"
        
    # Look for common fatty acid patterns
    fatty_chain = Chem.MolFromSmarts('C[CH2][CH2][CH2]')
    if not mol.HasSubstructMatch(fatty_chain):
        return False, "No typical fatty acid chain pattern found"

    # Success case
    features = []
    if alkene_count > 0:
        features.append(f"{alkene_count} double bonds")
    features.append(f"{epoxide_matches} epoxide ring(s)")
    if oxygen_count > 3:
        features.append(f"{oxygen_count-3} additional oxygen-containing groups")
    
    return True, f"Epoxy fatty acid with {carbon_count} carbons, " + ", ".join(features)