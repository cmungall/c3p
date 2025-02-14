"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamin
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, DataStructs

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a B vitamin based on its SMILES string.
    The B vitamins include vitamin B1 (thiamine), B2 (riboflavin), B3 (niacin),
    B5 (pantothenic acid), B6 (pyridoxine, pyridoxal, pyridoxamine),
    B7 (biotin), B9 (folic acid), and B12 (cobalamin).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a B vitamin, False otherwise
        str: Reason for classification
    """
    
    # List of known B vitamin SMILES strings
    b_vitamins = {
        'Vitamin B1 (Thiamine)': 'CC1=C(C)C=CN=C1C[C@H](O)CO',
        'Vitamin B2 (Riboflavin)': 'CC1=C(C)C=C2N(C[C@H](O)[C@H](O)[C@H](O)CO)C3=NC(=O)NC(=O)C3=NC2=C1',
        'Vitamin B3 (Niacin)': 'OC(=O)C1=CN=CC=C1',
        'Vitamin B5 (Pantothenic acid)': 'CC(C)(CO)C(=O)NCCC(=O)O',
        'Vitamin B6 (Pyridoxine)': 'NC1=NC=C(C(=C1O)CO)C',
        'Vitamin B7 (Biotin)': 'O=C1NC(=O)N2[C@@](CS1)([H])CCCCC2',
        'Vitamin B9 (Folic acid)': 'NC1=NC2=C(N1)C(=O)NC(=O)C2=CC=C3C=CC=CC3N',
        # Vitamin B12 (Cobalamin) is complex; we check for cobalt atom
    }
    
    # Convert input SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Generate fingerprint for input molecule
    input_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    
    # Threshold for Tanimoto similarity
    similarity_threshold = 0.8
    
    # Compare input molecule with known B vitamins
    for vitamin_name, vitamin_smiles in b_vitamins.items():
        # Convert vitamin SMILES to molecule
        vitamin_mol = Chem.MolFromSmiles(vitamin_smiles)
        if vitamin_mol is None:
            continue  # Skip invalid vitamin SMILES
        
        # Generate fingerprint for vitamin molecule
        vitamin_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(vitamin_mol, radius=2, nBits=2048)
        
        # Calculate Tanimoto similarity
        similarity = DataStructs.TanimotoSimilarity(input_fp, vitamin_fp)
        
        # If similarity is above threshold, classify as B vitamin
        if similarity >= similarity_threshold:
            return True, f"Molecule matches {vitamin_name}"
    
    # Check for the presence of cobalt for vitamin B12
    contains_cobalt = any(atom.GetAtomicNum() == 27 for atom in mol.GetAtoms())
    if contains_cobalt:
        return True, "Molecule contains cobalt characteristic of Vitamin B12 (Cobalamin)"
    
    return False, "Molecule does not match any known B vitamin"