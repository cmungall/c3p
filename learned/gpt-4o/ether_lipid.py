"""
Classifies: CHEBI:64611 ether lipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid is similar to a glycerolipid but with one or more ether linkages
    at the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flexibly recognize glycerol backbone with ether linkage
    # Pattern must account for chirality, branching, and variability
    glycerol_ether_pattern = Chem.MolFromSmarts("COC[C@H]([O*])COP([O-])(=O)OCC[N+]([CH3])([CH3])[CH3]")  # Example pattern including phosphocholine
    if not mol.HasSubstructMatch(glycerol_ether_pattern):
        return False, "No ether-linked glycerol backbone found"
    
    # Must have at least one ether linkage (-O-C-) substituting an ester
    ether_linkage_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    ether_linkages = mol.GetSubstructMatches(ether_linkage_pattern)
    if len(ether_linkages) < 1:
        return False, "No sufficient ether linkage found substituting ester linkage"

    return True, "Contains ether linkage with glycerol backbone"

# Example usage
example_smiles = "P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)C=C)COCCCCCCCCCCCCCCCCCC)([O-])=O"
result, reason = is_ether_lipid(example_smiles)
print(result, reason)