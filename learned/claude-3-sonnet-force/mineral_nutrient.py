"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: CHEBI:33621 mineral nutrient
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    A mineral nutrient is an inorganic compound containing essential mineral elements
    that must be ingested and absorbed in adequate amounts for proper metabolic and
    structural functions in the human body.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral nutrient, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the compound is inorganic
    if any(bond.GetBondType() == Chem.BondType.SINGLE for bond in mol.GetBonds() if bond.GetBeginAtom().GetSymbol() == "C" and bond.GetEndAtom().GetSymbol() == "C"):
        return False, "Organic compound (contains C-C bonds)"

    if any(atom.GetSymbol() == "N" and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 3 for atom in mol.GetAtoms()):
        return False, "Organic compound (contains amines)"

    if any(atom.GetSymbol() == "O" and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 2 for atom in mol.GetAtoms()):
        return False, "Organic compound (contains alcohols or ethers)"

    # List of mineral nutrient elements
    mineral_nutrient_elements = ["Na", "K", "Ca", "Mg", "Fe", "Zn", "Cu", "Mn", "I", "Se", "Cr", "Mo", "F", "Cl", "P", "S", "Ba", "Al", "Si", "Ni", "Co", "Sn", "V", "Li", "La", "Cs", "Sb"]

    # Check for presence of mineral nutrient elements
    has_mineral_element = any(atom.GetSymbol() in mineral_nutrient_elements for atom in mol.GetAtoms())
    if not has_mineral_element:
        return False, "No mineral nutrient elements found"

    # Check molecular weight range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 50 or mol_wt > 1000:
        return False, "Molecular weight outside typical range for mineral nutrients"

    # Check for common structural patterns and functional groups
    mineral_nutrient_patterns = [
        Chem.MolFromSmarts("[Na,K,Ca,Mg,Ba,Al,Cs,Fe,Zn,Cu,Mn,I,Se,Cr,Mo,Ni,Co,Sn,V,Li,La,Sb]"),  # Cation
        Chem.MolFromSmarts("[O-,F-,Cl-,Br-,I-,N-,S-,P-]"),  # Anion
        Chem.MolFromSmarts("[O-][P+]([O-])([O-])=O"),  # Phosphate
        Chem.MolFromSmarts("[O-][S+]([O-])(=O)=O"),  # Sulfate
        Chem.MolFromSmarts("[O-][N+]([O-])=O"),  # Nitrate
        Chem.MolFromSmarts("[O-][C+](=O)[O-]"),  # Carbonate
        Chem.MolFromSmarts("[O-][C+](=O)"),  # Carboxylate
        Chem.MolFromSmarts("[Cl,Br,I]~[La,Cs,Sb]"),  # Metal halides
        Chem.MolFromSmarts("[O-]~[Si]~[O-]"),  # Silicon oxides
    ]

    has_mineral_pattern = any(mol.HasSubstructMatch(pattern) for pattern in mineral_nutrient_patterns)
    if not has_mineral_pattern:
        return False, "No common mineral nutrient structural patterns found"

    # Check charge balance
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != 0:
        return False, "Unbalanced charge"

    # Handle exceptions and edge cases
    # ...

    return True, "Contains mineral nutrient elements and structural patterns"