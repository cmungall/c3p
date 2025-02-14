"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    A saturated fatty acid is a monocarboxylic acid with a hydrocarbon chain containing no carbon-carbon multiple bonds.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for exactly one carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1]')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Found {len(carboxylic_acid_matches)} carboxylic acid groups, need exactly 1"
    
    # Exclude molecules containing atoms other than C, H, O
    allowed_atomic_numbers = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_numbers:
            return False, f"Contains atom other than C, H, O: {atom.GetSymbol()}"

    # Exclude molecules with amide bonds (C(=O)-N)
    amide_pattern = Chem.MolFromSmarts('C(=O)N')
    if mol.HasSubstructMatch(amide_pattern):
        return False, "Contains amide bond, indicating a peptide"

    # Check for any carbon-carbon double or triple bonds
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
            if bond.GetBondType() in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
                return False, "Contains carbon-carbon multiple bonds"

    # Count the number of carbon atoms (excluding the carboxylic acid carbon)
    carboxylic_acid_carbons = {match[0] for match in carboxylic_acid_matches}
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIdx() not in carboxylic_acid_carbons]
    num_carbons = len(carbon_atoms)
    
    # No minimum carbon count to include short-chain fatty acids like propionic acid
    
    # Saturated fatty acids may contain rings and branches, so no need to exclude them
    
    return True, "Molecule is a saturated fatty acid"