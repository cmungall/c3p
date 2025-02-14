"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
"""
Classifies: CHEBI:27310 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    An N-acyl-L-alpha-amino acid is an L-alpha-amino acid carrying an N-acyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for alpha-amino acid backbone (C-C-C with amino and carboxyl groups)
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[C@@H](N)([C@@H](C(=O)O))[C@@H]")
    if not mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return False, "No L-alpha-amino acid backbone found"
    
    # Look for N-acyl group (N-C(=O)-)
    n_acyl_pattern = Chem.MolFromSmarts("N[C@@H](C(=O))")
    if not mol.HasSubstructMatch(n_acyl_pattern):
        return False, "No N-acyl group found"
    
    # Check for stereochemistry (L-configuration)
    mol = AllChem.AssignAtomChiralTagsFromStructure(mol)
    chiral_centers = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED]
    if len(chiral_centers) < 2:
        return False, "Not enough chiral centers to determine L-configuration"
    
    chiral_atoms = [mol.GetAtomWithIdx(idx) for idx in chiral_centers]
    for atom in chiral_atoms:
        if atom.GetHybridization() == Chem.HybridizationType.SP3 and atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
            return False, "Stereochemistry not consistent with L-configuration"
    
    return True, "Contains L-alpha-amino acid backbone with N-acyl substituent"