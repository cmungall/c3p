"""
Classifies: CHEBI:28863 flavanones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    Flavanones have a 2-phenyl-3,4-dihydro-2H-1-benzopyran-4-one skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for the flavanone core structure
    # chromane ring with carbonyl on position 4 and saturation at positions 2 and 3
    # This is a tricky SMARTS to write as it needs to make sure that the ring is exactly a 6 member ring and a specific arrangement of atoms
    # The carbonyl is explicitly placed by [C]=O, and the O is placed within the 6-member ring by [O]
    # The aromatic ring is represented by 'c'
    # The 2,3-dihydro part by [CH2][CH]
    flavanone_core_smarts = "[cH]1[cH][cH][c]2[O][CH]([CH2][CH]2[C]=O)[c]1"
    flavanone_core_pattern = Chem.MolFromSmarts(flavanone_core_smarts)

    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "Flavanone core structure not found"
    
    # Check for phenyl at position 2 (next to the oxygen, which is bonded to the [CH2])
    phenyl_smarts = "[cH]1[cH][cH][cH][cH][cH]1"
    phenyl_pattern = Chem.MolFromSmarts(phenyl_smarts)
    
    matches = mol.GetSubstructMatches(flavanone_core_pattern)
    
    #get the first match of the substructure
    if matches:
        match = matches[0]

        #find the C atom which is bonded to the 2-position of the chromane
        for atom_index in match:
            atom = mol.GetAtomWithIdx(atom_index)
            if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 1:  #identifies carbon in the CH group (at the 2-position)
                
                for neighbor in atom.GetNeighbors():
                  if neighbor.GetSymbol() == 'C' and neighbor.GetTotalNumHs() == 0: # Identifies carbon from the C attached to phenyl at position 2
                    
                    attached_phenyl = False
                    for phenyl_neighbor in neighbor.GetNeighbors():
                      if phenyl_neighbor.GetSymbol() == 'C' and phenyl_neighbor.IsInRing() and phenyl_neighbor.GetTotalNumHs() == 1:
                          attached_phenyl = True
                          break
                    if not attached_phenyl:
                        return False, "Phenyl substituent not found at correct position"

                    
                    phenyl_mol = Chem.MolFromSmiles(smiles)

                    
                    if not AllChem.FindSubstructMatch(phenyl_mol,phenyl_pattern):
                      return False, "Phenyl substituent not found"


                    # Check for at least one OH or OCH3 on the phenyl group
                    hydroxylated_phenyl_smarts = "[cH]1[cH][cH][cH]([O])[cH][cH]1"
                    methoxy_phenyl_smarts = "[cH]1[cH][cH][cH]([O][CH3])[cH][cH]1"
                    hydroxylated_phenyl_pattern = Chem.MolFromSmarts(hydroxylated_phenyl_smarts)
                    methoxy_phenyl_pattern = Chem.MolFromSmarts(methoxy_phenyl_smarts)
                    
                    if not mol.HasSubstructMatch(hydroxylated_phenyl_pattern) and not mol.HasSubstructMatch(methoxy_phenyl_pattern):
                      return False, "Phenyl group is missing at least one hydroxyl or methoxy substituent"
                    
                    
                    break

                break

    return True, "Molecule is a flavanone"