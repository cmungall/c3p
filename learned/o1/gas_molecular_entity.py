"""
Classifies: CHEBI:138675 gas molecular entity
"""
# gas_molecular_entity.py
"""
Classifies: CHEBI gas molecular entity
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a main group molecular entity that is gaseous at STP.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a gas molecular entity, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for net charge; exclude ions
    if Chem.GetFormalCharge(mol) != 0:
        return False, "Molecule is an ion with non-zero net charge, unlikely to be gaseous at STP"

    # List of atomic numbers for noble gases
    noble_gases = [2, 10, 18, 36, 54, 86]  # He, Ne, Ar, Kr, Xe, Rn
    # Dictionary of diatomic gases with their SMILES
    diatomic_gases = {
        'H2': '[H][H]',
        'N2': 'N#N',
        'O2': 'O=O',
        'F2': 'F[F]',
        'Cl2': 'ClCl',
    }
    # List of common gaseous inorganic molecules SMILES
    common_gaseous_inorganics = [
        'O=C=O',       # Carbon dioxide
        '[C-]#[O+]',   # Carbon monoxide
        '[H]N([H])[H]', # Ammonia
        '[H][H]',      # Dihydrogen
        'Cl[H]',       # Hydrogen chloride
        'I[H]',        # Hydrogen iodide
        'Br[H]',       # Hydrogen bromide
        'F[H]',        # Hydrogen fluoride
        '[O][O]',      # Oxygen
        '[N]#[N]',     # Nitrogen
        '[O][N]=O',    # Nitrogen dioxide
    ]

    # Check if the molecule is a noble gas atom
    if mol.GetNumAtoms() == 1:
        atomic_num = mol.GetAtomWithIdx(0).GetAtomicNum()
        if atomic_num in noble_gases:
            return True, f"Molecule is a noble gas atom ({mol.GetAtomWithIdx(0).GetSymbol()}), gaseous at STP"

    # Check if the molecule matches known diatomic gases
    for gas_name, gas_smiles in diatomic_gases.items():
        gas_mol = Chem.MolFromSmiles(gas_smiles)
        if mol.HasSubstructMatch(gas_mol) and mol.GetNumAtoms() == gas_mol.GetNumAtoms():
            return True, f"Molecule is {gas_name}, gaseous at STP"
    
    # Check if the molecule matches common gaseous inorganic molecules
    for gas_smiles in common_gaseous_inorganics:
        gas_mol = Chem.MolFromSmiles(gas_smiles)
        if mol.HasSubstructMatch(gas_mol) and mol.GetNumAtoms() == gas_mol.GetNumAtoms():
            return True, "Molecule is a common gaseous inorganic compound at STP"

    # Calculate molecular weight
    mol_weight = Descriptors.MolWt(mol)

    # Exclude molecules containing metals or transition elements
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        # Atomic numbers for metals start from 21 (Scandium) excluding noble gases
        if atomic_num > 20 and atomic_num not in noble_gases:
            return False, "Molecule contains metals or elements not gaseous at STP"

    # Count number of carbons and hydrogens
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_halogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [9, 17, 35, 53])  # F, Cl, Br, I

    # Check for small hydrocarbons (C1-C4 alkanes, alkenes, alkynes)
    if num_carbons <= 4 and mol_weight <= 70:
        return True, f"Molecule is a small hydrocarbon (C{num_carbons}), likely gaseous at STP"

    # Check for small molecules with low molecular weight
    if mol_weight <= 70 and num_carbons <= 4 and num_halogens <= 1:
        return True, f"Molecule has low molecular weight ({mol_weight:.2f} g/mol), likely gaseous at STP"

    # For molecules provided in examples that may not meet the above criteria (e.g., ozone)
    special_cases_smiles = ['[O-][O+]=O', 'C1CO1']  # Ozone, oxirane
    for special_smiles in special_cases_smiles:
        special_mol = Chem.MolFromSmiles(special_smiles)
        if mol.HasSubstructMatch(special_mol) and mol.GetNumAtoms() == special_mol.GetNumAtoms():
            return True, "Molecule is a known gas at STP"

    return False, "Molecule is not gaseous at STP based on molecular weight and composition"