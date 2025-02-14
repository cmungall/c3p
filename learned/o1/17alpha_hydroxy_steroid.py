"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
"""
Classifies: 17alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    A 17alpha-hydroxy steroid is a steroid with an alpha-oriented hydroxyl group at the C17 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure molecule has stereochemistry information
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    # Steroid backbone SMARTS with labeled atoms (C17 is labeled as atom 17)
    steroid_smarts = """
    [#6]-1=[#6]-[#6]-[#6]-2-[#6]-[#6]-3=[#6]-[#6]-[#6]-[#6]-[#6]-3-
    [#6]-[#6]-2-[#6]-1
    """

    # Create a template steroid molecule with labeled atoms
    steroid_template = Chem.MolFromSmiles('C1CCC2C3CCC4=CC(=O)CC[C@]4(C)[C@@H]3CC[C@]12C')
    if steroid_template is None:
        return False, "Failed to create steroid template"

    # Label the C17 atom in the template (atom index 16 in zero-based indexing)
    atom_map = {atom.GetIdx(): atom.GetIdx() + 1 for atom in steroid_template.GetAtoms()}
    steroid_template.GetAtomWithIdx(16).SetAtomMapNum(17)  # C17

    # Perform substructure match to align the molecule to the steroid template
    match = mol.GetSubstructMatch(steroid_template)
    if not match:
        return False, "Molecule does not match steroid backbone"

    # Create an atom map from template to molecule
    atom_map = {template_idx: mol_idx for template_idx, mol_idx in enumerate(match)}

    # Get the C17 atom in the molecule
    c17_idx = atom_map.get(16)  # Template atom index 16 corresponds to C17
    if c17_idx is None:
        return False, "C17 atom not found in molecule"

    c17_atom = mol.GetAtomWithIdx(c17_idx)

    # Check for hydroxyl group at C17
    has_oh = False
    for neighbor in c17_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 8:
            # Check if the oxygen is part of a hydroxyl group
            if len(neighbor.GetNeighbors()) == 1:
                has_oh = True
                oxygen_idx = neighbor.GetIdx()
                break

    if not has_oh:
        return False, "No hydroxyl group attached to C17"

    # Check stereochemistry at C17 (alpha orientation)
    stereo = c17_atom.GetChiralTag()
    if stereo == Chem.CHI_UNSPECIFIED:
        return False, "Chirality at C17 is unspecified"

    # Determine if the configuration is alpha at C17
    # For steroids, alpha orientation corresponds to S configuration at C17
    # RDKit uses CIP rules to assign R/S configuration
    c17_symbol = c17_atom.GetSymbol()
    if c17_atom.HasProp('_CIPCode'):
        cip_code = c17_atom.GetProp('_CIPCode')
        if cip_code != 'S':
            return False, f"Hydroxyl group at C17 is not alpha-oriented (CIP code: {cip_code})"
    else:
        return False, "Unable to assign CIP code to C17 atom"

    return True, "Molecule is a 17alpha-hydroxy steroid with correct stereochemistry at C17"