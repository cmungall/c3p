"""
Classifies: CHEBI:26167 polar amino acid
"""
"""
Classifies: CHEBI:26976 polar amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

# Define SMARTS patterns for amino acid backbone and polar groups
L_AA_BACKBONE = Chem.MolFromSmarts("[NH3+][C@@H](C(=O)[O-])C")  # L-amino acid backbone (zwitterionic)
D_AA_BACKBONE = Chem.MolFromSmarts("[NH3+][C@H](C(=O)[O-])C")  # D-amino acid backbone (zwitterionic)
POLAR_GROUPS = [
    Chem.MolFromSmarts("[OH1]"),  # Hydroxyl group
    Chem.MolFromSmarts("[NH2]"),  # Primary amine
    Chem.MolFromSmarts("[NH1]"),  # Secondary amine
    Chem.MolFromSmarts("[SH1]"),  # Thiol group
    Chem.MolFromSmarts("[NH+]=[NH2+]"),  # Guanidino group
    Chem.MolFromSmarts("c1cnc[nH]1"),  # Imidazole ring (non-aromatic)
    Chem.MolFromSmarts("c1nccn1"),  # Imidazole ring (aromatic)
]

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    A polar amino acid is an amino acid whose side chain is capable of forming
    one or more hydrogen bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polar amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for amino acid backbone (zwitterionic forms)
    if mol.HasSubstructMatch(L_AA_BACKBONE) or mol.HasSubstructMatch(D_AA_BACKBONE):
        # Remove the amino acid backbone
        side_chain = Chem.DeleteSubstructs(mol, L_AA_BACKBONE)
        side_chain = Chem.DeleteSubstructs(side_chain, D_AA_BACKBONE)

        # Check for polar groups in the side chain
        has_polar_group = any(side_chain.HasSubstructMatch(pattern) for pattern in POLAR_GROUPS)

        if has_polar_group:
            return True, "Contains amino acid backbone and polar side chain"
        else:
            return False, "Non-polar side chain"
    else:
        return False, "No amino acid backbone found"