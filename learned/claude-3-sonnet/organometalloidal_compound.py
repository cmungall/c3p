"""
Classifies: CHEBI:143084 organometalloidal compound
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    An organometalloidal compound has bonds between metalloid atoms (primarily As, but also B, Si, Ge, Sb, Te)
    and carbon atoms of an organyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organometalloidal compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define metalloids (focusing primarily on As since it's most common in examples)
    metalloids = [33]  # As atomic number
    
    # Track if we found valid metalloid-carbon bonds
    found_metalloid_carbon = False
    
    # Examine each atom in the molecule
    for atom in mol.GetAtoms():
        # Check if atom is a metalloid
        if atom.GetAtomicNum() in metalloids:
            # Look for carbon neighbors
            carbon_neighbors = [n for n in atom.GetNeighbors() 
                              if n.GetAtomicNum() == 6]
            
            if carbon_neighbors:
                # Check if any carbon is part of an organyl group
                for carbon in carbon_neighbors:
                    # Skip if carbon is part of CO3 or similar inorganic groups
                    oxygen_count = sum(1 for n in carbon.GetNeighbors() 
                                     if n.GetAtomicNum() == 8)
                    if oxygen_count >= 3:
                        continue
                        
                    # Valid metalloid-carbon bond found
                    found_metalloid_carbon = True
                    break
    
    if not found_metalloid_carbon:
        return False, "No valid metalloid-carbon bonds found"
        
    # Additional checks for arsenic compounds
    as_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 33]
    if not as_atoms:
        return False, "No arsenic atoms found"
        
    # Count total number of As-C bonds
    total_as_c_bonds = sum(len([n for n in as_atom.GetNeighbors() 
                               if n.GetAtomicNum() == 6]) 
                          for as_atom in as_atoms)
    
    # Verify that we have at least one As-C bond
    if total_as_c_bonds == 0:
        return False, "No As-C bonds found"
        
    # Success case
    return True, f"Found {total_as_c_bonds} As-C bond(s) in organometalloidal compound"