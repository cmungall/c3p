"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: CHEBI:38336 cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    A cannabinoid is a diverse group of pharmacologically active secondary metabolite 
    characteristic to Cannabis plant as well as produced naturally in the body by humans and animals. 
    Cannabinoids contain oxygen as a part of the heterocyclic ring or in the form of various functional groups. 
    They are subdivided on the basis of their origin.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for pyran, benzopyran, or indole ring systems
    pyran_pattern = Chem.MolFromSmarts("O1CCCCC1")
    benzopyran_pattern = Chem.MolFromSmarts("c1ccoc1")
    indole_pattern = Chem.MolFromSmarts("c1c[nH]cc1")
    if not (mol.HasSubstructMatch(pyran_pattern) or
            mol.HasSubstructMatch(benzopyran_pattern) or
            mol.HasSubstructMatch(indole_pattern)):
        return False, "No pyran, benzopyran, or indole ring found"

    # Check for long aliphatic chains (>5 carbons)
    aliphatic_pattern = Chem.MolFromSmarts("CCCCCC")
    if not mol.HasSubstructMatch(aliphatic_pattern):
        return False, "No long aliphatic chain found"

    # Check for oxygen-containing functional groups
    oxygen_pattern = Chem.MolFromSmarts("[#8]")
    if not mol.HasSubstructMatch(oxygen_pattern):
        return False, "No oxygen-containing functional group found"

    # Check for stereochemistry
    if mol.GetNumAtoms() > 10 and not any(atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for atom in mol.GetAtoms()):
        return False, "Molecule contains stereochemistry but it is not specified"

    return True, "Contains pyran/benzopyran/indole ring, long aliphatic chain, and oxygen-containing functional group"