"""
Classifies: CHEBI:25903 peptide antibiotic
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a peptide antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if molecule is organic
    if not mol.GetAtomWithIdx(0).GetSymbol() == "C":
        return False, "Not an organic molecule"
        
    # Calculate molecular properties
    mw = Descriptors.MolWt(mol)
    n_rotatable = Descriptors.NumRotatableBonds(mol)
    h_bond_donor = Lipinski.NumHDonors(mol)
    h_bond_acceptor = Lipinski.NumHAcceptors(mol)
    
    # Check if molecule meets basic peptide criteria
    if mw < 500 or mw > 5000:
        return False, f"Molecular weight ({mw:.2f} Da) out of typical peptide range"
    if n_rotatable < 5:
        return False, "Too few rotatable bonds for a peptide"
    if h_bond_donor < 2 or h_bond_acceptor < 4:
        return False, "Insufficient hydrogen bond donors/acceptors for a peptide"
        
    # Look for peptide bonds
    peptide_pattern = Chem.MolFromSmarts("C(=O)NCC")
    if not mol.HasSubstructMatch(peptide_pattern):
        return False, "No peptide bonds found"
    
    # Look for amino acid residues
    amino_acid_patterns = [Chem.MolFromSmarts(s) for s in [
        "[NH2]CCC(=O)O", # Basic
        "N[C@H](C)C(=O)O", # Ala
        "NC(C(=O)O)CC", # Asp/Asn
        # Add patterns for other amino acids here
    ]]
    
    amino_acid_matches = []
    for pattern in amino_acid_patterns:
        amino_acid_matches.extend(mol.GetSubstructMatches(pattern))
    
    if not amino_acid_matches:
        return False, "No amino acid residues found"
    
    # Look for non-ribosomal signature (ring systems)
    ring_info = mol.GetRingInfo()
    if not ring_info.AtomRings():
        return False, "No ring systems found (likely ribosomal peptide)"
    
    # Passed all checks, classify as peptide antibiotic
    return True, "Molecule meets criteria for peptide antibiotic"