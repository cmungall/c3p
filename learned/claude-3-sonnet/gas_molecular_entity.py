"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: CHEBI:36963 gas molecular entity
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

# Whitelist of known gas molecular entities at STP
GAS_WHITELIST = ['[220Rn]', '[219Rn]', '[222Rn]', '[He]', '[6He]', '[4He]', '[3He]', '[Kr]', '[Xe]', '[Ar]', '[Rn]', '[Ne]', '[C]']

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a gas molecular entity at standard temperature and pressure (STP; 0°C and 100 kPa).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gas molecular entity at STP, False otherwise
        str: Reason for classification
    """
    
    # Check whitelist
    if smiles in GAS_WHITELIST:
        return True, "Known gas molecular entity at STP"
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude charged species
    if any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms()):
        return False, "Contains charged atoms"
    
    # Exclude transition metals
    if any(atom.GetAtomicNum() in range(21, 31) or atom.GetAtomicNum() in range(39, 49) or atom.GetAtomicNum() in range(71, 81) for atom in mol.GetAtoms()):
        return False, "Contains transition metal atoms"
    
    # Check for common gas molecular entities
    if mol.GetNumAtoms() == 2:  # Diatomic gases
        atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
        if sorted(atoms) in [['Cl', 'Cl'], ['F', 'F'], ['H', 'H'], ['H', 'Cl'], ['H', 'I']]:
            return True, "Diatomic gas molecular entity"
    
    if mol.GetNumAtoms() == 3:  # Triatomic gases
        atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
        if sorted(atoms) == ['H', 'N', 'H'] or sorted(atoms) == ['O', 'O', 'O']:
            return True, "Triatomic gas molecular entity"
    
    # Check for small alkanes, alkenes, and alkynes
    if all(atom.GetAtomicNum() in [1, 6] for atom in mol.GetAtoms()):
        mol_wt = Descriptors.MolWt(mol)
        boiling_point = 1.96 * mol_wt ** 0.6 - 0.64  # Estimation based on Lee-Kesler method
        if boiling_point < 0:
            return True, "Small alkane, alkene, or alkyne gas molecular entity"
    
    return False, "Not a gas molecular entity at STP"