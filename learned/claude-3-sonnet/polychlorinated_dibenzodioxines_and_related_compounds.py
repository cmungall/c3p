"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: CHEBI:33631 polychlorinated dibenzodioxins and related compounds
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from typing import Tuple

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str) -> Tuple[bool, str]:
    """
    Determines if a molecule belongs to the class of polychlorinated dibenzodioxins and related compounds
    based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule belongs to the class, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of chlorine or bromine atoms
    has_halogens = any(atom.GetAtomicNum() in (17, 35) for atom in mol.GetAtoms())
    if not has_halogens:
        return False, "No chlorine or bromine atoms found"
    
    # Check for benzene rings
    benzene_rings = [ring for ring in mol.GetRingInfo().BondRings() if len(ring) == 6]
    if len(benzene_rings) < 2:
        return False, "Fewer than two benzene rings found"
    
    # Check for dioxin or furan rings
    dioxin_pattern = Chem.MolFromSmarts("O1c2ccccc2Oc3ccccc3O1")
    furan_pattern = Chem.MolFromSmarts("O1c2ccccc2oc3ccccc3O1")
    has_dioxin_furan = mol.HasSubstructMatch(dioxin_pattern) or mol.HasSubstructMatch(furan_pattern)
    
    # Check for biphenyl substructure
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1c2ccccc2")
    has_biphenyl = mol.HasSubstructMatch(biphenyl_pattern)
    
    # Check for halogenated benzene rings
    halogenated_benzene_pattern = Chem.MolFromSmarts("c1c(Cl)c(Cl)c(Cl)c(Cl)c1")
    has_halogenated_benzene = mol.HasSubstructMatch(halogenated_benzene_pattern)
    
    # Check for polybrominated substructures
    polybrominated_pattern = Chem.MolFromSmarts("c1c(Br)c(Br)c(Br)c(Br)c1")
    has_polybrominated = mol.HasSubstructMatch(polybrominated_pattern)
    
    if (has_dioxin_furan or has_biphenyl or has_halogenated_benzene or has_polybrominated) and has_halogens:
        return True, "Contains polychlorinated dibenzodioxin, dibenzofuran, halogenated biphenyl, or polybrominated substructures"
    else:
        return False, "Does not meet structural criteria for polychlorinated dibenzodioxins and related compounds"