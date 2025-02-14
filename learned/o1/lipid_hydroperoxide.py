"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
"""
Classifies: lipid hydroperoxide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is any lipid carrying one or more hydroperoxy (-OOH) substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to the molecule
    mol = Chem.AddHs(mol)

    # Check for hydroperoxy groups (-OOH)
    # Hydroperoxy group SMARTS pattern: [O][O][H]
    hydroperoxy_pattern = Chem.MolFromSmarts("[O][O][H]")
    hydroperoxy_matches = mol.GetSubstructMatches(hydroperoxy_pattern)
    num_hydroperoxy = len(hydroperoxy_matches)
    if num_hydroperoxy == 0:
        return False, "No hydroperoxy (-OOH) groups found"

    # Determine if the molecule is a lipid by checking for long hydrocarbon chains
    # Define a SMARTS pattern for a chain of at least 8 carbons
    long_chain_pattern = Chem.MolFromSmarts("[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long hydrocarbon chains found"

    # Optionally, check the total number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:
        return False, f"Too few carbon atoms ({c_count}) to be a lipid"

    return True, f"Contains {num_hydroperoxy} hydroperoxy group(s) and long hydrocarbon chains"

__metadata__ = {
    'chemical_class': {
        'name': 'lipid hydroperoxide',
        'definition': 'Any lipid carrying one or more hydroperoxy substituents.'
    }
}