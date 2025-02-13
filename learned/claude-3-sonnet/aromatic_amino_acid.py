"""
Classifies: CHEBI:33856 aromatic amino acid
"""
"""
Classifies: CHEBI:27561 aromatic amino acid
An amino acid whose structure includes an aromatic ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid is an amino acid that contains an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for amino acid functional group (-NH2 and -COOH)
    amino_pattern = Chem.MolFromSmarts("N")
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not amino_matches or not acid_matches:
        return False, "Missing amino acid functional groups"

    # Check for aromatic ring(s)
    aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() if Chem.Mol.Mol(mol).GetAromaticRingsAsMol().HasRingOnRing(ring)]
    if not aromatic_rings:
        return False, "No aromatic rings found"

    # Success
    return True, "Contains amino acid functional groups and aromatic ring(s)"