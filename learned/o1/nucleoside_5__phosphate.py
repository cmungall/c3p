"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate is a ribosyl or deoxyribosyl derivative of a pyrimidine
    or purine base in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns
    
    # Purine base pattern (adenine and guanine)
    purine_pattern = Chem.MolFromSmarts('n1c2ncnc2nc1')
    
    # Pyrimidine base pattern (cytosine, thymine, uracil)
    pyrimidine_pattern = Chem.MolFromSmarts('c1c[nH]c(=O)[nH]c1=O')
    
    # Ribose sugar pattern (five-membered ring with oxygen and hydroxyl groups)
    ribose_pattern = Chem.MolFromSmarts('[C@H]1(O)[C@@H](O)[C@H](O)[C@@H](O)O1')
    
    # Deoxyribose sugar pattern (ribose without hydroxyl at 2' position)
    deoxyribose_pattern = Chem.MolFromSmarts('[C@H]1(O)[C@@H](O)[C@H](O)[C@@H](H)O1')
    
    # Phosphate group pattern attached to primary carbon (5' position)
    phosphate_pattern = Chem.MolFromSmarts('[C;H2]-O-P(=O)(O)-O')
    
    # Check for purine or pyrimidine base
    has_purine = mol.HasSubstructMatch(purine_pattern)
    has_pyrimidine = mol.HasSubstructMatch(pyrimidine_pattern)
    
    if not (has_purine or has_pyrimidine):
        return False, "No purine or pyrimidine base found"
    
    # Check for ribose or deoxyribose sugar
    has_ribose = mol.HasSubstructMatch(ribose_pattern)
    has_deoxyribose = mol.HasSubstructMatch(deoxyribose_pattern)
    
    if not (has_ribose or has_deoxyribose):
        return False, "No ribose or deoxyribose sugar found"
    
    # Check for nucleotide linkage (base attached to sugar)
    # For simplicity, check if base and sugar are connected
    base_atoms = purine_pattern.GetSubstructMatch(mol) or pyrimidine_pattern.GetSubstructMatch(mol)
    sugar_atoms = ribose_pattern.GetSubstructMatch(mol) or deoxyribose_pattern.GetSubstructMatch(mol)
    
    if not base_atoms or not sugar_atoms:
        return False, "Base or sugar substructure not found"
    
    # Check if base is connected to sugar via N-glycosidic bond
    connected = False
    for atom_idx in base_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() in sugar_atoms:
                connected = True
                break
        if connected:
            break
    
    if not connected:
        return False, "Base is not connected to sugar via N-glycosidic bond"
    
    # Check for phosphate group attached to 5' carbon of sugar
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group attached to 5' carbon of sugar"
    
    # At this point, criteria are satisfied
    return True, "Molecule is a nucleoside 5'-phosphate"