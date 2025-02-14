"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
"""
Classifies: nucleoside phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate is a nucleobase-containing molecule where one or more of
    the sugar hydroxy groups has been converted into a mono- or poly-phosphate.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a nucleoside phosphate, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define common nucleobase patterns
    adenine_smarts = "c1nc2c(n1)nc(nc2N)N"
    guanine_smarts = "c1nc2c(n1)[nH]c(nc2O)N"
    cytosine_smarts = "C1=CN=CN=C1O"
    thymine_smarts = "C1=CN=CN=C1O"
    uracil_smarts = "O=C1NC=CC(=O)N1"
    nucleobase_patterns = [
        Chem.MolFromSmarts(adenine_smarts),
        Chem.MolFromSmarts(guanine_smarts),
        Chem.MolFromSmarts("C1=CN=CN=C1N"),  # Cytosine simplified pattern
        Chem.MolFromSmarts("C1=CN=CN=C1N"),  # Thymine simplified pattern
        Chem.MolFromSmarts(uracil_smarts)
    ]

    # Check for presence of nucleobase
    nucleobase_found = False
    for base_pattern in nucleobase_patterns:
        if mol.HasSubstructMatch(base_pattern):
            nucleobase_found = True
            break
    if not nucleobase_found:
        return False, "No nucleobase found"

    # Define sugar moiety pattern (ribose or deoxyribose)
    sugar_pattern = Chem.MolFromSmarts("C1C(C(C(O1)CO)O)O")  # Ribose ring with hydroxyl groups
    deoxysugar_pattern = Chem.MolFromSmarts("C1C(C(C(O1)CO)O)")  # Deoxyribose lacks one OH
    if not mol.HasSubstructMatch(sugar_pattern) and not mol.HasSubstructMatch(deoxysugar_pattern):
        return False, "No ribose or deoxyribose sugar moiety found"

    # Check for glycosidic bond between nucleobase and sugar
    glycosidic_bond_pattern = Chem.MolFromSmarts("c1ncnc1[C@H]2O[C@@H](C[C@@H](O)[C@H]2O)O")  # Simplified pattern
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond between nucleobase and sugar found"

    # Define phosphate group pattern attached to sugar's hydroxyl group
    phosphate_ester_pattern = Chem.MolFromSmarts("O[P](=O)(O)O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")  # Phosphate ester linked to sugar
    phosphate_patterns = [
        Chem.MolFromSmarts("OP(=O)(O)O"),  # Monophosphate
        Chem.MolFromSmarts("OP(=O)(OP(=O)(O)O)O"),  # Diphosphate
        Chem.MolFromSmarts("OP(=O)(OP(=O)(OP(=O)(O)O)O)O")  # Triphosphate
    ]

    # Check for phosphate groups attached to sugar
    phosphate_found = False
    for phosphate_pattern in phosphate_patterns:
        if mol.HasSubstructMatch(phosphate_pattern):
            phosphate_found = True
            break
    if not phosphate_found:
        return False, "No phosphate groups attached to the sugar"

    return True, "Molecule is a nucleoside phosphate with nucleobase, sugar, and phosphate group(s)"