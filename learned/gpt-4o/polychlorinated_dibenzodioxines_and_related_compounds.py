"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Classifies a molecule as a polychlorinated dibenzodioxin or related compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule belongs to the class, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for polychlorinated dibenzodioxins
    dioxin_pattern = Chem.MolFromSmarts('O1c2ccccc2Oc3ccccc13')
    
    # Define SMARTS patterns for dibenzofurans
    furan_pattern = Chem.MolFromSmarts('O1c2ccccc2-c3c1ccccc3')

    # Define SMARTS patterns for halogenated biphenyls (at least tetrachlorinated)
    biphenyl_pattern = Chem.MolFromSmarts('c1ccccc1-c1ccccc1')
    
    # Counting halogen atoms directly
    num_cl = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'Cl')
    num_br = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'Br')
    
    # Determine the presence of the key structure motifs
    has_dioxin = mol.HasSubstructMatch(dioxin_pattern)
    has_furan = mol.HasSubstructMatch(furan_pattern)
    has_biphenyl = mol.HasSubstructMatch(biphenyl_pattern)
    
    # Check whether the halogen count is significant
    high_halogens = num_cl + num_br >= 4  # At least tetrachlorinated compounds

    # Ensure at least one significant structure type and significant halogenation
    if (has_dioxin or has_furan or has_biphenyl) and high_halogens:
        return True, "Classified as polychlorinated dibenzodioxin or related compound"
    
    return False, "Does not meet classification criteria for polychlorinated dibenzodioxins or related compounds"

# Example usage
example_smiles = "Clc1ccc(cc1)-c1cc(Cl)c(Cl)cc1Cl"
result, reason = is_polychlorinated_dibenzodioxines_and_related_compounds(example_smiles)
print(result, reason)