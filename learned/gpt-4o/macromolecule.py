"""
Classifies: CHEBI:33839 macromolecule
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
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight too low ({mol_wt} < 1000 g/mol)"
    
    # Look for repeating units which could indicate a polymeric structure
    # Simplified approach: check for duplicate substructs within reasonable size
    fragment_counts = Chem.rdMolDescriptors.CalcNumRings(mol)
    if fragment_counts < 5:
        return False, "No significant repeating structural motifs"
    
    # Check for high complexity: high number of atoms and bonds
    num_atoms = mol.GetNumAtoms()
    num_bonds = mol.GetNumBonds()
    if num_atoms > 100 and num_bonds > 150:
        return True, f"Contains {num_atoms} atoms and {num_bonds} bonds indicating high complexity"
    
    # If both high molecular weight and repeating motifs are detected
    if mol_wt >= 1000 and fragment_counts >= 5:
        return True, f"Considerable molecular weight and repeating structural motifs"
    
    return False, "Structure does not meet macromolecule criteria"