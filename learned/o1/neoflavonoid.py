"""
Classifies: CHEBI:71971 neoflavonoid
"""
"""
Classifies: neoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is any 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for 1-benzopyran core with atom mapping
    benzopyran_smarts = "[c:1]1[c:2][c:3][c:4]2[o:5][c:6][c:7][c:8][c:9]2[c:10]1"
    benzopyran_mol = Chem.MolFromSmarts(benzopyran_smarts)
    matches = mol.GetSubstructMatches(benzopyran_mol)

    if not matches:
        return False, "No 1-benzopyran core detected"

    # Iterate over matches to find if any has an aryl substituent at position 4
    for match in matches:
        # Get the atom index of position 4
        atom_idx_pos4 = match[3]  # Atom mapping number 4 corresponds to index 3
        atom_pos4 = mol.GetAtomWithIdx(atom_idx_pos4)

        # Examine neighbors of position 4 atom that are not part of the benzopyran core
        for neighbor in atom_pos4.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in match:
                # Identify the substituent group attached at position 4
                # Check if it is an aryl group (aromatic ring connected directly)
                bond = mol.GetBondBetweenAtoms(atom_idx_pos4, neighbor_idx)
                if bond is None:
                    continue

                # Use atom index to get the fragment (substituent)
                substituent = Chem.PathToSubmol(mol, Chem.FindAtomEnvironmentOfRadiusN(mol, 3, neighbor_idx))

                # Define SMARTS pattern for aryl group (any aromatic ring)
                aryl_pattern = Chem.MolFromSmarts("a1aaaaa1")  # Six-membered aromatic ring
                if substituent.HasSubstructMatch(aryl_pattern):
                    return True, "Molecule is a neoflavonoid with aryl substituent at position 4"

                # Also check for five-membered and seven-membered aromatic rings
                aryl_pattern5 = Chem.MolFromSmarts("a1aaaa1")
                aryl_pattern7 = Chem.MolFromSmarts("a1aaaaaa1")
                if substituent.HasSubstructMatch(aryl_pattern5) or substituent.HasSubstructMatch(aryl_pattern7):
                    return True, "Molecule is a neoflavonoid with aryl substituent at position 4"

    return False, "No aryl substituent detected at position 4 of the benzopyran core"