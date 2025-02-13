"""
Classifies: CHEBI:26267 proanthocyanidin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    A proanthocyanidin is a flavonoid oligomer with two or more units of hydroxyflavans linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proanthocyanidin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for hydroxyflavan pattern (Common in all proanthocyanidins)
    hydroxyflavan_pattern = Chem.MolFromSmarts("c1cc(O)ccc1C2(CO)Oc3cc(O)ccc3C2")  # Simplified
    if not mol.HasSubstructMatch(hydroxyflavan_pattern):
        return False, "No hydroxyflavan units found"

    # Check for linkage pattern (4->8 or 4->6 linkages)
    linkage_pattern = Chem.MolFromSmarts("C1(COC2)CC(C3(OCC4))") # Simplified 4->8 linkage
    linkage_matches = mol.GetSubstructMatches(linkage_pattern)
    if len(linkage_matches) < 1:
        return False, "Missing flavan linkages (4->8 or 4->6)"

    # Check for oligomeric structure (multiple flavan units)
    flavan_units_pattern = Chem.MolFromSmarts("C1(CO)C2(OCC3)")  # Base flavan structure
    flavan_units_matches = mol.GetSubstructMatches(flavan_units_pattern)
    if len(flavan_units_matches) < 2:
        return False, f"Only {len(flavan_units_matches)} flavan unit(s) found, need at least 2"

    # Optional: Check for common modifications like gallate esters
    # gallate_pattern = Chem.MolFromSmarts("OC(=O)C1=CC(O)=C(O)C=C1") 
    # if any(mol.HasSubstructMatch(patt) for patt in [gallate_pattern]):
    #     print("Gallate ester groups detected")

    return True, "Contains multiple hydroxyflavan units with appropriate linkages"