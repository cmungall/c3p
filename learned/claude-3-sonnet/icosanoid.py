"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: CHEBI:38108 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    Icosanoids are signaling molecules derived from oxidation of C20 essential fatty acids
    like icosapentaenoic acid (EPA), arachidonic acid (AA), and dihomo-gamma-linolenic acid (DGLA).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define substructure patterns for icosanoid precursors
    epa_pattern = Chem.MolFromSmarts("CCCCCCCCC\C=C\C/C=C\C/C=C\CCCC(O)=O")  # EPA
    aa_pattern = Chem.MolFromSmarts("CCCCCCCCC\C=C\C/C=C\C/C=C\CCCC(O)=O")  # AA
    dgla_pattern = Chem.MolFromSmarts("CCCCCCCCC(O)\C=C\C/C=C\C/C=C\CCCC(O)=O")  # DGLA

    # Check for precursor substructures
    if mol.HasSubstructMatch(epa_pattern) or mol.HasSubstructMatch(aa_pattern) or mol.HasSubstructMatch(dgla_pattern):
        return True, "Contains a substructure derived from EPA, AA, or DGLA"

    # Define common icosanoid functional group patterns
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    epoxy_pattern = Chem.MolFromSmarts("[O;X2]1[C;X4][C;X4]1")
    keto_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[!#8]")
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][!#8]")

    # Check for common functional groups
    has_hydroxy = mol.HasSubstructMatch(hydroxy_pattern)
    has_epoxy = mol.HasSubstructMatch(epoxy_pattern)
    has_keto = mol.HasSubstructMatch(keto_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)

    if has_hydroxy or has_epoxy or has_keto or has_ester:
        return True, "Contains common icosanoid functional groups"

    return False, "No evidence of being an icosanoid"