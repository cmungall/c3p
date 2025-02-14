"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
"""
Classifies: CHEBI:36973 glycosaminoglycan
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    A glycosaminoglycan is any polysaccharide containing a substantial proportion
    of aminomonosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for patterns of aminomonosaccharide residues
    amino_sugar_pattern = Chem.MolFromSmarts("[N;X3;H2,H1][-;!@C]")
    amino_sugar_matches = mol.GetSubstructMatches(amino_sugar_pattern)

    # Check if there are multiple aminomonosaccharide residues
    if len(amino_sugar_matches) < 2:
        return False, "Insufficient aminomonosaccharide residues"

    # Check for polysaccharide backbone
    polysaccharide_pattern = Chem.MolFromSmarts("[OX2]~[CX4]~[OX2]")
    polysaccharide_matches = mol.GetSubstructMatches(polysaccharide_pattern)
    if len(polysaccharide_matches) < 2:
        return False, "No polysaccharide backbone found"

    # Check for additional glycosidic bonds
    glycosidic_pattern = Chem.MolFromSmarts("[OX2]~[CX4]~[OX2]~[CX4]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if not glycosidic_matches:
        return False, "No glycosidic bonds found"

    # Count nitrogen atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 2:
        return False, "Too few nitrogen atoms for glycosaminoglycan"

    # Check molecular weight - glycosaminoglycans typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for glycosaminoglycan"

    return True, "Contains aminomonosaccharide residues and a polysaccharide backbone"