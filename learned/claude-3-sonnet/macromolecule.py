"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Classifies: CHEBI:24116 macromolecule
A macromolecule is a molecule of high relative molecular mass, the structure of which essentially
comprises the multiple repetition of units derived, actually or conceptually, from molecules of low
relative molecular mass.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macromolecule, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight - macromolecules typically >1000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, "Molecular weight too low for macromolecule"
    
    # Look for repeating substructures (e.g. polymers)
    smarts_patterns = ["[$(C-N)]~[$(N-C)]", "[$(C-O)]~[$(O-C)]"]  # peptide and polysaccharide patterns
    for pattern in smarts_patterns:
        patt = Chem.MolFromSmarts(pattern)
        matches = mol.GetSubstructMatches(patt)
        if len(matches) > 2:  # more than 2 repeating units
            return True, "Contains repeating substructures indicative of a macromolecule"
    
    # If high MW but no repeating units, could be small biomolecule
    return False, "High molecular weight but no repeating substructures found"