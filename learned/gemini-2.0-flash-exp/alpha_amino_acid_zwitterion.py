"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid zwitterion based on its SMILES string.
    An alpha-amino acid zwitterion has a protonated amino group and a deprotonated carboxyl group on the alpha carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid zwitterion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the basic substructure of zwitterion
    zwitterion_pattern = Chem.MolFromSmarts("[CH]([NH3+])C([O-])=O")
    if not mol.HasSubstructMatch(zwitterion_pattern):
        return False, "Missing zwitterionic alpha-amino acid core structure"

    # Verify that there is at least one [NH3+] and one [O-]
    nh3_pattern = Chem.MolFromSmarts("[NH3+]")
    o_minus_pattern = Chem.MolFromSmarts("[O-]")
    if not mol.HasSubstructMatch(nh3_pattern) or not mol.HasSubstructMatch(o_minus_pattern):
        return False, "Missing either positively charged nitrogen or negatively charged oxygen"


    # Exclude peptides (more than one alpha aminoacid)
    zwitterion_count = len(mol.GetSubstructMatches(zwitterion_pattern))
    if zwitterion_count > 1:
        return False, "Molecule contains multiple alpha-amino acid units (possibly a peptide)"


    return True, "Molecule contains the zwitterionic alpha-amino acid structure"