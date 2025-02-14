"""
Classifies: CHEBI:33856 aromatic amino acid
"""
"""
Classifies: CHEBI:33567 aromatic amino acid
An amino acid whose structure includes an aromatic ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.

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
    
    # Check for amino acid pattern
    amino_acid_pattern = Chem.MolFromSmarts("[NX3H2][CX4H]([CX3](=O)[OX1-,OX2H1])[CX4H]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid pattern found"
    
    # Check for aromatic rings
    aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() if mol.GetRingAtoms(ring).IsAromatic()]
    if not aromatic_rings:
        return False, "No aromatic rings found"
    
    # Check if aromatic ring is part of the molecule backbone
    aromatic_ring_atoms = set.union(*[set(mol.GetRingAtoms(ring)) for ring in aromatic_rings])
    amino_acid_atoms = set(mol.GetSubstructMatches(amino_acid_pattern)[0])
    if not aromatic_ring_atoms.intersection(amino_acid_atoms):
        return False, "Aromatic ring not part of backbone"
    
    return True, "Contains an aromatic ring and an amino acid backbone"