"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid (MIA) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is an MIA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for indole substructure: benzene ring fused to a pyrrole (indole)
    indole_pattern = Chem.MolFromSmarts('c1c[cH]cc2c1[nH]cc2')
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "Indole moiety not found"

    # Look for characteristic terpenoid features (such as cyclized isoprene groups)
    terpene_pattern = Chem.MolFromSmarts('C(C(=O)O)C(C)C')  # Pseudo pattern for terpene linkage
    if not mol.HasSubstructMatch(terpene_pattern):
        return False, "Terpene-like linkage not found"

    # Examine complexity based on molecular weight and presence of oxygen
    is_large = Descriptors.MolWt(mol) > 350  # MIAs are generally large
    has_ester_oxygen = any(atom.GetAtomicNum() == 8 and atom.GetDegree() == 2 for atom in mol.GetAtoms())
    
    # Chiral centers indicating complex natural products
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    
    if not (is_large and has_ester_oxygen and len(chiral_centers) >= 3):
        return False, "Does not meet expected structure complexity for MIAs"

    return True, "Contains indole moiety with complex MIA-specific features and ester linkage"

# Example usage:
smiles = "C/C=C\1/CN2CC[C@]34C5=CC(C=6C=C7C(=CC6OC)[C@]89CCN%10C/C(=C\C)/[C@](CCC8%10N7C)([C@]9(C(=O)OC)[H])[H])=C(C=C5N(C)[C@]4(C2C[C@@]1([C@@]3(C(=O)OC)[H])[H])[H])OC"
print(is_monoterpenoid_indole_alkaloid(smiles))