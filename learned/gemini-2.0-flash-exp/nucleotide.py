"""
Classifies: CHEBI:36976 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide is a nucleoside phosphate resulting from the condensation of the 3 or 5
    hydroxy group of a nucleoside with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one phosphate group
    phosphate_pattern = Chem.MolFromSmarts("[P](=[OX1])([OX2])([OX2])")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for a pentose sugar - flexible pattern 5 C or 4C and one O
    pentose_pattern1 = Chem.MolFromSmarts("[C;X4]~[C;X4]~[C;X4]~[C;X4]~[C;X4]") # 5 C
    pentose_pattern2 = Chem.MolFromSmarts("[C;X4]~[C;X4]~[C;X4]~[O;X2]~[C;X4]") # 4 C and one O
    if not (mol.HasSubstructMatch(pentose_pattern1) or mol.HasSubstructMatch(pentose_pattern2)):
       return False, "No pentose sugar ring found"

    # Check for a nucleobase - an N atom in an aromatic cycle
    nucleobase_pattern = Chem.MolFromSmarts("[n;r]")
    if not mol.HasSubstructMatch(nucleobase_pattern):
        return False, "No nucleobase found"
    
    # Check if a phosphate is connected to an oxygen.
    phosphate_connection_pattern = Chem.MolFromSmarts("[P](=[OX1])([OX2])([OX2])~[OX2]")
    if not mol.HasSubstructMatch(phosphate_connection_pattern):
        return False, "Phosphate not connected to an oxygen"
    
    # Check that pentose sugar is connected to an oxygen
    pentose_sugar_connection_pattern = Chem.MolFromSmarts("[#6]~[OX2]")
    if not mol.HasSubstructMatch(pentose_sugar_connection_pattern):
        return False, "No hydroxyl group found"

    # Check for N-glycosidic bond
    N_glycosidic_pattern = Chem.MolFromSmarts("[N;!H0]~[#6]~[O;X2]")
    if not mol.HasSubstructMatch(N_glycosidic_pattern):
       return False, "No N-glycosidic linkage"
    
    return True, "Contains a nucleobase, a pentose sugar and a phosphate group connected correctly."