"""
Classifies: CHEBI:16389 ubiquinones
"""
from rdkit import Chem

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    Ubiquinones typically have a 2,3-dialkoxy benzoquinone core and a polyisoprenyl side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 2,3-dialkoxy benzoquinone core
    benzoquinone_core_pattern = Chem.MolFromSmarts("COc1c(OC)cc(=O)cc1=O")
    if not mol.HasSubstructMatch(benzoquinone_core_pattern):
        return False, "No 2,3-dialkoxy benzoquinone core found"

    # Look for long isoprenoid chain
    isoprenoid_pattern = Chem.MolFromSmarts("C=C(C)C")
    isoprenoid_matches = mol.GetSubstructMatches(isoprenoid_pattern)
    if len(isoprenoid_matches) < 3:  # Typically ubiquinones have 6-10 units
        return False, "Insufficient isoprenoid units, found " + str(len(isoprenoid_matches))

    return True, "Contains 2,3-dialkoxy benzoquinone core with sufficient isoprenoid units"

# Test cases for debugging
test_smiles = [
    "COC1=C(OC)C(=O)C(C\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)=C(C)C1=O",  # ubiquinone-9
    "O=C1C(OC)=C(OC)C(=O)C=C1C(C(O)C)C",  # 5-(3-hydroxybutan-2-yl)-2,3-dimethoxycyclohexa-2,5-diene-1,4-dione
]

for smiles in test_smiles:
    result, reason = is_ubiquinones(smiles)
    print(f"SMILES: {smiles}\nIs Ubiquinone: {result}\nReason: {reason}\n")