"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import EditableMol

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon with one or more hydrogen atoms replaced by a nitro group (-NO2).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a nitrohydrocarbon, False otherwise.
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for invalid atoms. Only C, H, N, O are allowed
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in [1, 6, 7, 8]:  # 1(H), 6(C), 7(N), 8(O)
            return False, f"Contains invalid atom: {atom.GetSymbol()}"

    # 2. Check for at least one nitro group
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    if not mol.HasSubstructMatch(nitro_pattern):
        return False, "No nitro group found"

    # 3. Create a copy of the molecule, remove the nitro groups and check that what remains is a hydrocarbon
    mol_copy = Chem.Mol(mol)  # create a copy of the molecule
    edit_mol = EditableMol(mol_copy) # create an editable molecule

    nitro_matches = mol.GetSubstructMatches(nitro_pattern)

    if not nitro_matches:
         return False, "No nitro group found"
    
    atoms_to_remove_idxs = [] # collect the indices of the atoms to remove, they are shared by the original and the copy
    bonds_to_remove = []
    
    for match in nitro_matches:
        # get all the atoms involved in the substructure match
        atoms_to_remove_idx = []
        for idx in match:
            atoms_to_remove_idx.append(mol_copy.GetAtomWithIdx(idx).GetIdx())
            
        nitro_nitrogen_idx = -1
        for atom_idx in match:
          if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 7:
            nitro_nitrogen_idx = atom_idx
            break
        
        if nitro_nitrogen_idx == -1:
          return False, "Could not find nitrogen in -NO2 group"

        carbon_neighbor_count = 0;
        carbon_neighbor_idx = -1;
        for neighbor in mol.GetAtomWithIdx(nitro_nitrogen_idx).GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                carbon_neighbor_count+=1
                carbon_neighbor_idx = neighbor.GetIdx()

        if carbon_neighbor_count != 1:
             return False, "Nitro group is not bonded to exactly one carbon"
        
        atoms_to_remove_idxs.extend(atoms_to_remove_idx)
        
        # remove bonds inside of the nitro group
        for bond in mol_copy.GetBonds():
            if bond.GetBeginAtomIdx() in atoms_to_remove_idx and bond.GetEndAtomIdx() in atoms_to_remove_idx:
              bonds_to_remove.append(bond.GetIdx())

    bonds_to_remove = sorted(list(set(bonds_to_remove)), reverse=True)

    edit_mol.BeginBatchEdit()
    for bond_idx in bonds_to_remove:
        try:
            edit_mol.RemoveBond(mol_copy.GetBondWithIdx(bond_idx).GetBeginAtomIdx(), mol_copy.GetBondWithIdx(bond_idx).GetEndAtomIdx())
        except:
            return False, "Could not remove bond from molecule copy"
    
    atoms_to_remove_idxs = sorted(list(set(atoms_to_remove_idxs)), reverse=True)
    for atom_idx in atoms_to_remove_idxs:
        try:
             edit_mol.RemoveAtom(atom_idx)
        except:
            return False, "Could not remove atom from molecule copy"
    edit_mol.CommitBatchEdit()

    mol_copy = edit_mol.GetMol()


    for atom in mol_copy.GetAtoms():
      if atom.GetAtomicNum() not in [1,6]:
        return False, "After removing nitro groups, non-hydrocarbon fragment remains"

    # 4. Check that all nitrogens are part of a -NO2 group
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen
            if atom.GetTotalDegree() != 3: # must be bound to 3 other atoms
                return False, "Nitrogen is not part of -NO2 group (does not have 3 bonds)"
            
            oxygen_count = 0;
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    oxygen_count += 1
            if oxygen_count != 2: # must be bonded to two oxygens
                return False, "Nitrogen is not part of -NO2 group (does not have 2 oxygen neighbors)"

    return True, "Molecule is a nitrohydrocarbon"