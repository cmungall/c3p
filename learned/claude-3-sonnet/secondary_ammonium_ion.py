"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: CHEBI:36836 secondary ammonium ion
An organic cation obtained by protonation of any secondary amino compound; major species at pH 7.3.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for charged atoms
    charged_atoms = [atom for atom in mol.GetAtoms() if atom.GetFormalCharge() != 0]
    if not charged_atoms:
        return False, "No charged atoms found"

    # Find the charged nitrogen atom
    charged_n = None
    for atom in charged_atoms:
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
            charged_n = atom
            break
    if charged_n is None:
        return False, "No charged nitrogen atom found"

    # Check if the charged nitrogen is attached to two carbon atoms and a hydrogen
    neighbors = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in charged_n.GetNeighbors()]
    if neighbors.count('C') != 2 or 'H' not in neighbors:
        return False, "Charged nitrogen does not meet secondary amine criteria"

    # Check if the molecule contains only C, H, N, O, S, P, F, Cl, Br, I
    allowed_atoms = [6, 1, 7, 8, 16, 15, 9, 17, 35, 53]  # C, H, N, O, S, P, F, Cl, Br, I
    if any(atom.GetAtomicNum() not in allowed_atoms for atom in mol.GetAtoms()):
        return False, "Molecule contains disallowed atoms"

    return True, "Molecule contains a protonated secondary amine group"