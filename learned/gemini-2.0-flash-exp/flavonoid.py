"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid has a 1-benzopyran core with an aryl substituent at position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for the 1-benzopyran core, allowing a double bond in the pyran ring
    benzopyran_core_smarts_1 = "[c]1[c][c][c]([o][c]2[c][c][c][c][c]2)[c]1" # 1-benzopyran with the o in position 1 (relative to the aromatic ring), single bond in pyran
    benzopyran_core_smarts_2 = "[c]1[c][c][c]([o,C]=[c]2[c][c][c][c][c]2)[c]1" # 1-benzopyran with the o in position 1 (relative to the aromatic ring), double bond in pyran
    
    # Define SMARTS for the aryl substituent at position 2. The aryl group is directly linked to the [c] next to the oxygen.
    aryl_substituent_smarts = "[c]1[c][c][c][c][c]1" # phenyl ring

    # Convert SMARTS to Mol objects
    benzopyran_core_mol_1 = Chem.MolFromSmarts(benzopyran_core_smarts_1)
    benzopyran_core_mol_2 = Chem.MolFromSmarts(benzopyran_core_smarts_2)
    aryl_substituent_mol = Chem.MolFromSmarts(aryl_substituent_smarts)

    # Check for the benzopyran core
    if not (mol.HasSubstructMatch(benzopyran_core_mol_1) or mol.HasSubstructMatch(benzopyran_core_mol_2)):
        return False, "No 1-benzopyran core found"
    
    # Find the core atoms (to determine the position 2)
    core_match_1 = mol.GetSubstructMatch(benzopyran_core_mol_1)
    core_match_2 = mol.GetSubstructMatch(benzopyran_core_mol_2)

    if core_match_1:
        # Position of the aryl group is the atom next to the O (that is not the aromatic ring atom)
         for match_index in mol.GetSubstructMatches(benzopyran_core_mol_1):
            
            o_index = match_index[3]
            # find the carbon neighbor to oxygen that is not aromatic
            for neighbor in mol.GetAtomWithIdx(o_index).GetNeighbors():
                if neighbor.GetIdx() in match_index:
                    continue
                
                c_index = neighbor.GetIdx()
                
                # Check for the aryl group
                
                for match_aryl_index in mol.GetSubstructMatches(aryl_substituent_mol):

                    for aryl_atom in match_aryl_index:
                        if mol.GetBondBetweenAtoms(c_index, aryl_atom):
                            return True, "Has 1-benzopyran core and an aryl substituent at position 2"

            return False, "No aryl substituent at position 2"

    elif core_match_2:
         for match_index in mol.GetSubstructMatches(benzopyran_core_mol_2):
            
            o_index = match_index[3] # position of oxygen
            # find the carbon neighbor to oxygen that is not aromatic
            for neighbor in mol.GetAtomWithIdx(o_index).GetNeighbors():
                if neighbor.GetIdx() in match_index:
                    continue
                
                c_index = neighbor.GetIdx()
                
                # Check for the aryl group
                
                for match_aryl_index in mol.GetSubstructMatches(aryl_substituent_mol):

                    for aryl_atom in match_aryl_index:
                        if mol.GetBondBetweenAtoms(c_index, aryl_atom):
                            return True, "Has 1-benzopyran core and an aryl substituent at position 2"
            return False, "No aryl substituent at position 2"

    return False, "Should not reach here"