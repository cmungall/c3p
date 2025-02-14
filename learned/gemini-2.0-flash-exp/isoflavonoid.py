"""
Classifies: CHEBI:50753 isoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is a 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Define the 1-benzopyran core with a SMARTS pattern, including the possibility for carbonyl and other substitutions at position 4
    # The wildcard [*:1] is used for later substructure query
    benzopyran_pattern = Chem.MolFromSmarts("c1ccccc2[o,C]([*:1])[c](=[O,C])[c]21")
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "Not a 1-benzopyran"

    # Get the indices of the matching atoms, the position that needs to be substituted is index 1 in the SMARTS pattern
    match = mol.GetSubstructMatch(benzopyran_pattern)
    if not match:
        return False, "No match for benzopyran core"
    c3_index = match[1]

    # 3. Check for an aryl substituent (phenyl or substituted phenyl) at C3
    #   Here we use a generic aryl as a substituent, to allow for the variations of the substituent
    aryl_pattern = Chem.MolFromSmarts("[cX3]1[cX3][cX3][cX3][cX3][cX3]1")
    
    # Get neighbors for c3 atom
    c3_atom = mol.GetAtomWithIdx(c3_index)
    neighbors = [neighbor.GetIdx() for neighbor in c3_atom.GetNeighbors()]

    # Now we check that at least one of these neighbors matches the aryl_pattern
    is_isoflavonoid = False
    reason = "Substituent at position 3 is not an aryl"
    for neighbor_idx in neighbors:
        sub_mol = Chem.RWMol(mol)
        # Keep only the atoms in the neighborhood, to not have too many false positives
        atom_to_keep = [c3_index, neighbor_idx] + [x.GetIdx() for x in mol.GetAtomWithIdx(neighbor_idx).GetNeighbors()]
        for i in reversed(range(sub_mol.GetNumAtoms())):
            if i not in atom_to_keep:
                sub_mol.RemoveAtom(i)
        sub_mol = Chem.Mol(sub_mol)
        if sub_mol.HasSubstructMatch(aryl_pattern):
            is_isoflavonoid = True
            reason = "1-Benzopyran with an aryl substituent at position 3"
            break

    return is_isoflavonoid, reason