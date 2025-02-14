"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
"""
Classifies: O-acyl-L-carnitine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    An O-acyl-L-carnitine is an O-acylcarnitine in which the carnitine component has L-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the quaternary ammonium group pattern
    quat_nitrogen_smarts = "[N+](C)(C)C"
    quat_nitrogen_mol = Chem.MolFromSmarts(quat_nitrogen_smarts)
    if not mol.HasSubstructMatch(quat_nitrogen_mol):
        return False, "No quaternary ammonium group found"

    # Define the carboxylate group pattern
    carboxylate_smarts = "C(=O)[O-]"
    carboxylate_mol = Chem.MolFromSmarts(carboxylate_smarts)
    if not mol.HasSubstructMatch(carboxylate_mol):
        return False, "No carboxylate group found"

    # Define the esterified hydroxy group at position 3 with correct stereochemistry
    # L-carnitine has (S) configuration, which corresponds to [C@@H] in SMILES/SMARTS
    ester_smarts = "[C@@H](COC(=O)[#6])[CH2][N+](C)(C)C"
    ester_mol = Chem.MolFromSmarts(ester_smarts)
    if not mol.HasSubstructMatch(ester_mol, useChirality=True):
        return False, "Esterified hydroxy group with correct stereochemistry not found"

    # Check if the acyl group is attached via ester bond
    ester_bond_found = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6) or (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 8):
                # Check for C-O single bond which could be part of ester
                neighb_atoms = [atom.GetAtomicNum() for atom in [atom1, atom2]]
                if 6 in neighb_atoms and 8 in neighb_atoms:
                    ester_bond_found = True
                    break
    if not ester_bond_found:
        return False, "No ester bond found at hydroxy group"

    return True, "Molecule is an O-acyl-L-carnitine with correct stereochemistry and esterified hydroxy group"

__metadata__ = {   
    'chemical_class': {   
        'id': None,
        'name': 'O-acyl-L-carnitine',
        'definition': 'An O-acylcarnitine in which the carnitine component has L-configuration.',
    },
    'config': {   'llm_model_name': None,
                  'f1_threshold': None,
                  'max_attempts': None,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': None,
                  'max_negative_in_prompt': None,
                  'max_instances_in_prompt': None,
                  'test_proportion': None},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
}