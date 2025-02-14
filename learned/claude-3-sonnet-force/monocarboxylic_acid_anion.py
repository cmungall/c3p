"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
"""
Classifies: CHEBI:36363 monocarboxylic acid anion

A monocarboxylic acid anion is formed when the carboxy group of a monocarboxylic acid is deprotonated. It should contain exactly one carboxylate (-COO-) group, and may have additional carbonyl and hydroxy groups. Zwitterionic structures should be excluded.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for zwitterionic structures (exclude positively charged atoms)
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() > 0:
            return False, "Molecule contains positively charged atoms (zwitterionic structure)"

    # Count carboxylate anions (-COO-) and carbonyl oxygens (C=O)
    carboxylate_count = 0
    carbonyl_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "O" and atom.GetFormalCharge() == -1:
            # Negatively charged oxygen (carboxylate anion)
            carboxylate_count += 1
        elif atom.GetSymbol() == "O" and atom.GetIsAromatic() == False:
            # Non-aromatic oxygen (carbonyl oxygen)
            carbonyl_count += 1

    # Criteria for monocarboxylic acid anion
    if carboxylate_count == 1 and carbonyl_count >= 1:
        return True, "Molecule contains one carboxylate anion and at least one carbonyl oxygen"
    elif carboxylate_count > 1:
        return False, f"Molecule contains {carboxylate_count} carboxylate anions (should be 1)"
    else:
        return False, "Molecule does not contain a carboxylate anion"