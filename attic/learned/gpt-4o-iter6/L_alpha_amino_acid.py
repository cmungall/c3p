"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False, useLegacyImplementation=False)

    # Define a more refined SMARTS pattern for an L-alpha-amino acid: central C with H, N, and C=O attached
    # and optically (stereo) active (R/S), check L configuration via @
    alpha_aa_smarts = '[C@@H](N)C(=O)O'
    alpha_aa_pattern = Chem.MolFromSmarts(alpha_aa_smarts)

    if not mol.HasSubstructMatch(alpha_aa_pattern):
        return False, "No L-alpha-amino acid backbone pattern found"

    # Check if there's at least one 'S' configuration at the chiral center for L-amino acids
    # Acknowledge that RDKit uses Cahn-Ingold-Prelog (CIP) priority with R/S configuration
    for center, chirality in chiral_centers:
        if chirality == 'R':
            atom = mol.GetAtomWithIdx(center)
            if atom.GetSymbol() == 'C':
                # Verify that the expected neighbors form valid L-alpha-amino acid
                n_neighbors = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
                if 'N' in n_neighbors and 'O' in n_neighbors:
                    return True, "Chiral center with L-configuration identified, amino and carboxyl groups confirmed."

    return False, "No suitable chiral configuration or functional groups for L-alpha-amino acid found."