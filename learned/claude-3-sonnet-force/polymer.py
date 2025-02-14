"""
Classifies: CHEBI:60027 polymer
"""
Here's a program to classify chemical entities as polymers based on their SMILES string:

"""
Classifies: CHEBI:36962 polymer

A polymer is a mixture, which is composed of macromolecules of different kinds 
and which may be differentiated by composition, length, degree of branching etc.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    A polymer is a large molecule composed of repeating structural units (monomers).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polymer, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for multiple connected components (salts, mixtures, etc.)
    components = Chem.GetMolFrags(mol)
    if len(components) > 1:
        # Polymers can be mixtures, so proceed with analysis
        pass
    
    # Look for repeating substructures
    patterns = Chem.FindMolChemicalEnvironmentOfAtoms(mol, AllChem.SSM.SSMatSmin())
    num_patterns = len(patterns)
    
    if num_patterns == 1:
        return False, "Only one unique structural unit (monomer) found"
    
    # Calculate molecular weight and count monomers
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    monomer_wts = [pattern.weight for pattern in patterns]
    monomer_counts = [int(mol_wt / wt) for wt in monomer_wts]
    
    # Polymers typically have high MW and repeating units
    if mol_wt < 500 or max(monomer_counts) < 5:
        return False, "Molecular weight too low or insufficient repeating units for polymer"
    
    # Check for long carbon chains (common in polymers)
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if carbon_chain_matches:
        return True, "Contains long carbon chains and repeating structural units (polymeric)"
    
    return True, "Contains repeating structural units (polymeric)"