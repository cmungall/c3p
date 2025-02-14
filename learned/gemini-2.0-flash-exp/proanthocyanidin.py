"""
Classifies: CHEBI:26267 proanthocyanidin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    A proanthocyanidin is a flavonoid oligomer formed by condensation of two or more hydroxyflavan units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proanthocyanidin, False otherwise.
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core flavan unit pattern
    flavan_pattern = Chem.MolFromSmarts("c1cc(O)c(C[C@H](O)[C@H]2Oc3cc(O)cc(O)c2c3)cc1O")
    if flavan_pattern is None:
        return False, "Invalid flavan unit SMARTS pattern"

    flavan_matches = mol.GetSubstructMatches(flavan_pattern)
    num_flavan_units = len(flavan_matches)
    if num_flavan_units < 2:
        return False, f"Found {num_flavan_units} flavan units, need at least 2."

    # Check for linkages between units. This is a complex problem, let's try to find two flavan units connected by single bond using a smarts pattern that contains a connecting bond and two flavan units
    
    # Define a pattern for two flavan units connected by a single bond
    linked_flavan_pattern = Chem.MolFromSmarts("[c1cc(O)c(C[C@H](O)[C@H]2Oc3cc(O)cc(O)c2c3)cc1O]~[c1cc(O)c(C[C@H](O)[C@H]2Oc3cc(O)cc(O)c2c3)cc1O]")
    if linked_flavan_pattern is None:
        return False, "Invalid linked flavan unit SMARTS pattern"

    linked_flavan_matches = mol.GetSubstructMatches(linked_flavan_pattern)

    # Check if there are two flavan units directly linked together
    if len(linked_flavan_matches) < 1 :
        # Look for patterns with C-O-C linkages between units.
        linked_flavan_pattern_O = Chem.MolFromSmarts("[c1cc(O)c(C[C@H](O)[C@H]2Oc3cc(O)cc(O)c2c3)cc1O]~O~[c1cc(O)c(C[C@H](O)[C@H]2Oc3cc(O)cc(O)c2c3)cc1O]")
        if linked_flavan_pattern_O is None:
            return False, "Invalid linked flavan unit SMARTS pattern"
        linked_flavan_matches_O = mol.GetSubstructMatches(linked_flavan_pattern_O)
        if len(linked_flavan_matches_O) < 1:
            return False, "Not linked flavan units"

    #Additional check: Count the number of aromatic rings. Proanthocyanidins should have a minimum of 4
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if aromatic_rings < 4:
        return False, "Too few aromatic rings for a proanthocyanidin"

    # Check the molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Too small molecular weight for a proanthocyanidin"


    return True, "Contains at least two linked flavan units."