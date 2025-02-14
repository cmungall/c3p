"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: CHEBI:37492 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide is a nucleoside phosphate resulting from the condensation of the 3 or 5 hydroxy
    group of a nucleoside with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for nucleoside backbone (purine or pyrimidine ring with sugar attached)
    nucleoside_pattern = Chem.MolFromSmarts("[N&R1]1[C&R2]=[N&R1][C&R2]=[N&R1][C&R2]=[N&R1]1[C&R2]*[C@H]2[C@H]([C@H]([C@H](O2)O)O)O")
    if not mol.HasSubstructMatch(nucleoside_pattern):
        return False, "No nucleoside backbone found"
    
    # Look for phosphate group (-OP(O)(O)=O)
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"
    
    # Check if phosphate is attached to sugar at 3' or 5' position
    sugar_atom_idx = nucleoside_pattern.GetAtomsMatchingSmarts("[C@H]2[C@H]([C@H]([C@H](O2)O)O)O")[0]
    sugar_atom = mol.GetAtomWithIdx(list(nucleoside_pattern.GetAtomsMatchingSmarts("[C@H]2[C@H]([C@H]([C@H](O2)O)O)O"))[0])
    phosphate_attached = False
    for bond in sugar_atom.GetBonds():
        if bond.GetOtherAtom(sugar_atom).GetSymbol() == "P":
            phosphate_attached = True
            break
    if not phosphate_attached:
        return False, "Phosphate not attached to sugar at 3' or 5' position"
    
    return True, "Contains nucleoside backbone with phosphate group attached"