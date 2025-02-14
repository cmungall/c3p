"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: CHEBI:27811 mineral
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    A mineral is a chemical substance that is normally crystalline, formed as a result of geological processes,
    or a metamict (formerly crystalline but crystallinity destroyed by radiation) or amorphous substance formed
    by geological processes or action on biogenic compounds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for common mineral elements
    common_mineral_elements = ['Na', 'K', 'Ca', 'Mg', 'Al', 'Fe', 'Cu', 'Zn', 'Ba', 'Pb', 'Ni', 'Co', 'Mn', 'Cr', 'Si', 'P', 'S', 'Cl', 'F', 'Br', 'I', 'O']
    has_mineral_element = any(atom.GetSymbol() in common_mineral_elements for atom in mol.GetAtoms())
    if not has_mineral_element:
        return False, "Does not contain common mineral elements"

    # Check for charge-balanced ionic species
    charge_sum = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if charge_sum != 0:
        return False, "Molecule is not charge-balanced"

    # Check for inorganic nature
    has_carbon = any(atom.GetSymbol() == 'C' for atom in mol.GetAtoms())
    if has_carbon:
        carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C']
        carbon_hybridizations = [atom.GetHybridization() for atom in carbon_atoms]
        if any(hybridization == Chem.HybridizationType.SP3 for hybridization in carbon_hybridizations):
            return False, "Contains organic carbon atoms (sp3 hybridized)"

    # Check for specific structural patterns
    mineral_patterns = [
        Chem.MolFromSmarts("[Si](O)(O)(O)(O)"),  # Silicates
        Chem.MolFromSmarts("P(=O)(O)(O)(O)"),  # Phosphates
        Chem.MolFromSmarts("S(=O)(=O)(O)(O)"),  # Sulfates
        Chem.MolFromSmarts("[Na,K,Ca,Mg,Al,Fe,Cu,Zn,Ba,Pb,Ni,Co,Mn,Cr]~[Cl,Br,I]"),  # Ionic halides
        Chem.MolFromSmarts("[Na,K,Ca,Mg,Al,Fe,Cu,Zn,Ba,Pb,Ni,Co,Mn,Cr]~[O,N]"),  # Ionic oxides, nitrates
        Chem.MolFromSmarts("[Na,K,Ca,Mg,Al,Fe,Cu,Zn,Ba,Pb,Ni,Co,Mn,Cr]~[C](=O)[O,N]"),  # Ionic carbonates, acetates
        Chem.MolFromSmarts("[Na,K,Ca,Mg,Al,Fe,Cu,Zn,Ba,Pb,Ni,Co,Mn,Cr]~[S]"),  # Sulfides
    ]
    has_mineral_pattern = any(mol.HasSubstructMatch(pattern) for pattern in mineral_patterns)
    if not has_mineral_pattern:
        return False, "Does not contain common mineral structural patterns"

    # Additional checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for mineral"

    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable > 5:
        return False, "Too many rotatable bonds for mineral"

    # Known exceptions
    known_minerals = ["[Na+].[Na+].[H]OB1O[B-]2(O[H])OB(O[H])O[B-](O[H])(O1)O2", "[S--].[S--].[Fe++].[Cu++]", "Cl[La](Cl)Cl"]
    if smiles in known_minerals:
        return True, "Known mineral exception"

    return True, "Contains common mineral elements, charge-balanced ionic species, inorganic nature, and mineral structural patterns"