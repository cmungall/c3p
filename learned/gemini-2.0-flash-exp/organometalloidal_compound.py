"""
Classifies: CHEBI:143084 organometalloidal compound
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    An organometalloidal compound has bonds between one or more metalloid atoms (B, Si, Ge, As, Sb, Te) 
    and one or more carbon atoms of an organyl group. The metalloid atom must be bound to an organic carbon.

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

    # Define metalloid atomic numbers
    metalloid_atoms = [5, 14, 32, 33, 51, 52]  # B, Si, Ge, As, Sb, Te

    # Check if any metalloid atoms exist in the molecule
    has_metalloid = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in metalloid_atoms:
            has_metalloid = True
            break
    if not has_metalloid:
        return False, "No metalloid atoms found"
    
    # Check for metalloid-carbon bonds where the carbon is also part of an organic group (sp3 hybridization check)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in metalloid_atoms:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Check if neighbor is carbon
                    
                    # Get carbon hybridization
                    hybridization = neighbor.GetHybridization()
                   
                    if hybridization == Chem.rdchem.HybridizationType.SP3:
                        return True, "Contains at least one metalloid-carbon bond in an organic group"

    return False, "No metalloid-carbon bond to an organic group found"