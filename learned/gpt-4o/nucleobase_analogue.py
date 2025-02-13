"""
Classifies: CHEBI:67142 nucleobase analogue
"""
from rdkit import Chem

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is defined as a molecule that can substitute for a normal nucleobase in nucleic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Patterns for nucleobase-like heterocycles
    pyrimidine_patterns = [
        Chem.MolFromSmarts("c1[nH]cnc1"),  # Pyrimidine with at least one imino group
        Chem.MolFromSmarts("c1nc[nH]c1"),  # Pyrimidine with another nitrogen configuration
    ]
    purine_patterns = [
        Chem.MolFromSmarts("c1ncnc2[nH]c[nH]c12"),  # Purine-like with two imino groups
    ]
    
    # Check for pyrimidine or purine-like structures
    for pyrimidine in pyrimidine_patterns:
        if mol.HasSubstructMatch(pyrimidine):
            # Ensure it has modifications typical for analogues
            halogen_pattern = Chem.MolFromSmarts("[F,Cl,Br,I]")
            hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
            if mol.HasSubstructMatch(halogen_pattern) or mol.HasSubstructMatch(hydroxyl_pattern):
                return True, "Contains pyrimidine-like heterocycle with modifications"
            
    for purine in purine_patterns:
        if mol.HasSubstructMatch(purine):
            return True, "Contains purine-like heterocycle"
        
    # Look for typical nucleobase-like functional group arrangements
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")  # Carbonyl group
    amino_pattern = Chem.MolFromSmarts("[NX3H2,X3H1]")  # Primary and secondary amines

    if mol.HasSubstructMatch(carbonyl_pattern) and mol.HasSubstructMatch(amino_pattern):
        # Additional check ensures the presence of at least a core heterocycle typical of nucleobases
        if any(mol.HasSubstructMatch(p) for p in pyrimidine_patterns + purine_patterns):
            return True, "Contains key nucleobase functional groups with heterocycle"
    
    return False, "Does not match features of known nucleobase analogues"