"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
"""
Classifies: CHEBI:35738 semisynthetic derivative
'Any organic molecular entity derived from a natural product by partial chemical synthesis.'
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFingerprintGenerator

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on its SMILES string.
    A semisynthetic derivative is any organic compound derived from a natural product
    by partial chemical synthesis.

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

    # Calculate fingerprint and look for natural product substructures
    fp = AllChem.GetMorganFingerprint(mol, 3)
    natural_product_bits = [
        fp.GetBitInfo(bit).PathFromAtom(mol.GetAtomWithIdx(0), atomCounts=mol.GetNumAtoms(), strict=False, atomMapNumbers=True) 
        for bit in fp.GetNonzeroPositions()
    ]
    natural_product_substructs = [Chem.MolFromSmarts(bit) for bit in natural_product_bits]

    # Check if the molecule contains any natural product substructures
    has_natural_product_substructs = any(mol.HasSubstructMatch(substruct) for substruct in natural_product_substructs)

    # Check if molecule is synthetic (not a natural product) by looking for synthetic building blocks
    synthetic_smarts = ['[#7-&!N1,N2,N3]', '[#6&a]', '[#16]', '[#9]', '[#17]', '[#35]', '[#53]']
    synthetic_substructs = [Chem.MolFromSmarts(smart) for smart in synthetic_smarts]
    is_synthetic = any(mol.HasSubstructMatch(substruct) for substruct in synthetic_substructs)

    if has_natural_product_substructs and is_synthetic:
        return True, "Contains natural product substructures and synthetic building blocks, likely a semisynthetic derivative"
    elif has_natural_product_substructs:
        return False, "Appears to be a natural product, no synthetic building blocks found"
    else:
        return False, "No natural product substructures found, likely fully synthetic"