"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: Epoxy Fatty Acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid is a heterocyclic fatty acid containing an epoxide ring as part of its structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for epoxide ring (three-membered cyclic ether)
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide ring found"

    # Optionally, check for long hydrocarbon chain typical of fatty acids
    # Count the number of carbon atoms excluding the carboxylic acid group
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(c_atoms)
    if c_count < 12:
        return False, f"Molecule has {c_count} carbon atoms, which is less than the typical fatty acid length"

    return True, "Molecule contains both a carboxylic acid group and an epoxide ring characteristic of epoxy fatty acids"