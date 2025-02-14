"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: ribonucleoside
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside is a nucleoside where the sugar component is D-ribose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure molecule has explicit hydrogens
    mol = Chem.AddHs(mol)

    # Define ribose sugar pattern (β-D-ribofuranose with correct stereochemistry)
    ribose_smarts = """
    [C@H]1(O)[O][C@@H]([C@H](O)[C@@H]1O)CO
    """
    ribose_pattern = Chem.MolFromSmarts(ribose_smarts)
    if ribose_pattern is None:
        return False, "Error in ribose pattern"

    # Find matches for ribose sugar
    ribose_matches = mol.GetSubstructMatches(ribose_pattern, useChirality=True)
    if not ribose_matches:
        return False, "No D-ribose sugar moiety found"

    # Define nucleobase patterns
    nucleobase_patterns = []

    # Adenine pattern
    adenine_smarts = 'n1c[nH]c2c1ncnc2N'
    adenine_pattern = Chem.MolFromSmarts(adenine_smarts)
    nucleobase_patterns.append(adenine_pattern)

    # Guanine pattern
    guanine_smarts = 'n1c(=O)[nH]c2c1ncnc2N'
    guanine_pattern = Chem.MolFromSmarts(guanine_smarts)
    nucleobase_patterns.append(guanine_pattern)

    # Cytosine pattern
    cytosine_smarts = 'n1c(N)nc(=O)cc1'
    cytosine_pattern = Chem.MolFromSmarts(cytosine_smarts)
    nucleobase_patterns.append(cytosine_pattern)

    # Uracil pattern
    uracil_smarts = 'O=C1NC(=O)C=C1'
    uracil_pattern = Chem.MolFromSmarts(uracil_smarts)
    nucleobase_patterns.append(uracil_pattern)

    # Hypoxanthine pattern (found in inosine)
    hypoxanthine_smarts = 'n1c(=O)[nH]c2c1ncnc2'
    hypoxanthine_pattern = Chem.MolFromSmarts(hypoxanthine_smarts)
    nucleobase_patterns.append(hypoxanthine_pattern)

    # Check for nucleobase attached to sugar
    nucleobase_found = False
    for pattern in nucleobase_patterns:
        if pattern is None:
            continue
        # Check for the nucleobase pattern in the molecule
        base_matches = mol.GetSubstructMatches(pattern)
        if not base_matches:
            continue
        # Check for connectivity between sugar and base (N-glycosidic bond)
        for base_match in base_matches:
            base_atoms = set(base_match)
            for ribose_match in ribose_matches:
                ribose_atoms = set(ribose_match)
                # Find bonds between ribose and base
                for atom_idx in ribose_atoms:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() in base_atoms:
                            # Check if the bond is between anomeric carbon (C1') and nitrogen of the base
                            if atom.GetAtomicNum() == 6 and nbr.GetAtomicNum() == 7:
                                nucleobase_found = True
                                break
                    if nucleobase_found:
                        break
                if nucleobase_found:
                    break
            if nucleobase_found:
                break
        if nucleobase_found:
            break

    if not nucleobase_found:
        return False, "No nucleobase moiety attached to D-ribose sugar via N-glycosidic bond"

    return True, "Contains D-ribose sugar moiety attached to nucleobase via N-glycosidic bond"