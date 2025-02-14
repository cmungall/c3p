"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
"""
Classifies: semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on its SMILES string.
    A semisynthetic derivative is an organic molecule derived from a natural product by partial chemical synthesis.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a semisynthetic derivative, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule is organic (contains carbon atoms)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule is not an organic compound"

    # Compute molecular weight
    mol_weight = Descriptors.ExactMolWt(mol)

    # Compute number of rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()

    # Compute number of chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    num_chiral_centers = len(chiral_centers)

    # Criteria for natural product core
    if num_rings < 3 or num_chiral_centers < 2 or mol_weight < 500:
        return False, "Molecule may not contain a complex natural product core"

    # Generate Murcko Scaffold (core structure of molecule)
    scaffold = Chem.Scaffolds.MurckoScaffold.GetScaffoldForMol(mol)
    scaffold_smiles = Chem.MolToSmiles(scaffold, isomericSmiles=True)

    # Check if scaffold is significantly smaller than molecule (indicating modifications)
    mol_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    if scaffold_smiles == mol_smiles:
        return False, "Molecule appears to be a natural product without modifications"

    # Check for synthetic modifications
    synthetic_modifications = False

    # Expanded list of synthetic functional groups uncommon in natural products
    synthetic_patterns = {
        'Halogenation': '[F,Cl,Br,I]',
        'Nitro group': '[NX3](=O)=O',
        'Azide': 'N=[N+]=[N-]',
        'Alkyne': 'C#C',
        'Trifluoromethyl': 'C(F)(F)F',
        'Sulfonamide': 'S(=O)(=O)N',
        'Sulfonic acid': 'S(=O)(=O)[O-]',
        'Protecting group': '[Si]',
        'Unnatural phosphate': 'P(=O)(O)(O)',
        'Sulfonate ester': 'S(=O)(=O)O',
        'Propargyl': 'C#CC',
        'Aryl ether': 'c-O-c',
    }

    for mod_name, smarts in synthetic_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            synthetic_modifications = True
            break  # One modification is enough to consider

    if not synthetic_modifications:
        return False, "No synthetic modifications detected"

    return True, "Molecule contains complex natural product core with evidence of synthetic modifications"