"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule is a polychlorinated dibenzodioxine, dibenzofuran, or polychlorinated/polybrominated biphenyl.
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if the molecule belongs to the class, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count halogens (Cl or Br)
    halogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [17, 35])
    if halogen_count < 2:
        return False, "Must have at least two chlorine or bromine atoms."

    # Check for dibenzodioxin core structure
    dioxin_pattern = Chem.MolFromSmarts("c1ccccc1-O-c1ccccc1-O-c1ccccc1")
    furan_pattern = Chem.MolFromSmarts("c1ccccc1-O-c1ccccc1-c1ccccc1")
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c1ccccc1")


    if mol.HasSubstructMatch(dioxin_pattern):
        return True, "Polychlorinated dibenzodioxin detected"
    if mol.HasSubstructMatch(furan_pattern):
          return True, "Polychlorinated dibenzofuran detected"
    if mol.HasSubstructMatch(biphenyl_pattern):
        return True, "Polychlorinated or polybrominated biphenyl detected"
    
    # Check for edge cases based on a few examples
    ambigol_pattern = Chem.MolFromSmarts("c1c(O)cc(Cl)c(O[c2c(Cl)cc(Cl)cc2])c1Cl")
    formicamycin_pattern = Chem.MolFromSmarts("C1=C(O)C(Cl)=C(O)C2=C1C([C@H]3CC=4C(Cl)=C(OC)C=C(C4C([C@]3(C2=O)O)=O)")
    if mol.HasSubstructMatch(ambigol_pattern):
      return True, "Ambigol-related structure detected."
    if mol.HasSubstructMatch(formicamycin_pattern):
      return True, "Formicamycin-related structure detected"

    return False, "Does not match any of the defined substructures or it lacks sufficient chlorines or bromines"