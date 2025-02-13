"""
Classifies: CHEBI:48030 tetrapeptide
"""
"""
Classifies: CHEBI:36357 tetrapeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide is defined as any molecule that contains four amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for 4 peptide bonds
    peptide_bond = Chem.MolFromSmarts("C(=O)NCC")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond)
    if len(peptide_bond_matches) != 4:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, expected 4"
    
    # Check for 4 amino acid residues
    aa_pattern = Chem.MolFromSmarts("N[C@H](C)C(=O)")  # Simple pattern for alpha-amino acids
    aa_matches = mol.GetSubstructMatches(aa_pattern)
    if len(aa_matches) != 4:
        return False, f"Found {len(aa_matches)} amino acid residues, expected 4"
    
    # Check molecular weight - tetrapeptides typically 300-800 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 800:
        return False, "Molecular weight outside typical range for tetrapeptide"
    
    # Check elemental composition - tetrapeptides should contain N, O, C, H and possibly S
    allowed_atoms = set([6, 7, 8, 16])  # C, N, O, S
    atoms = set([atom.GetAtomicNum() for atom in mol.GetAtoms()])
    if not atoms.issubset(allowed_atoms):
        return False, "Found unexpected atoms for tetrapeptide"
    
    return True, "Contains 4 amino acid residues connected by peptide bonds"