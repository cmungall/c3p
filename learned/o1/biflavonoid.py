"""
Classifies: CHEBI:50128 biflavonoid
"""
"""
Classifies: biflavonoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    A biflavonoid is a flavonoid dimer obtained by the coupling of at least two flavonoid units (aryl-substituted benzopyran rings or derivatives), resulting in two ring systems being joined together by a single atom or bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a biflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general flavonoid core SMARTS pattern (C6-C3-C6 skeleton)
    flavonoid_core_smarts = '[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1~[#6]2~[#6]~[#6]~[#6]~[#6]~[#8]~2'  # General flavonoid core
    flavonoid_core = Chem.MolFromSmarts(flavonoid_core_smarts)
    if flavonoid_core is None:
        return False, "Invalid flavonoid core SMARTS pattern"

    # Find matches for the flavonoid core
    flavonoid_matches = mol.GetSubstructMatches(flavonoid_core)
    num_flavonoid_cores = len(flavonoid_matches)

    if num_flavonoid_cores < 2:
        return False, f"Found {num_flavonoid_cores} flavonoid core(s), need at least 2"

    # Create sets of atoms for each flavonoid core
    flavonoid_core_atoms = [set(match) for match in flavonoid_matches]

    # Check if flavonoid cores are connected via a single atom or bond
    found_connected = False
    for i in range(num_flavonoid_cores):
        for j in range(i+1, num_flavonoid_cores):
            # Check for direct bonds or shared atoms between cores i and j
            inter_core_bonds = 0
            inter_core_atoms = flavonoid_core_atoms[i] & flavonoid_core_atoms[j]
            if inter_core_atoms:
                inter_core_bonds = len(inter_core_atoms)
            else:
                for atom_i in flavonoid_core_atoms[i]:
                    for atom_j in flavonoid_core_atoms[j]:
                        bond = mol.GetBondBetweenAtoms(atom_i, atom_j)
                        if bond is not None:
                            inter_core_bonds += 1
            if inter_core_bonds == 1:
                found_connected = True
                break
        if found_connected:
            break

    if not found_connected:
        return False, "Flavonoid cores are not connected via a single atom or bond"

    return True, "Contains at least two flavonoid cores connected via a single atom or bond"


__metadata__ = {   'chemical_class': {   'name': 'biflavonoid',
                              'definition': 'A flavonoid oligomer that is obtained by the oxidative coupling of at least two units of aryl-substituted benzopyran rings or its substituted derivatives, resulting in the two ring systems being joined together by a single atom or bond.'},
        'config': {   'llm_model_name': 'lbl/claude-sonnet'},
        'success': True}