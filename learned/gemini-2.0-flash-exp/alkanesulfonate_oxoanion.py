"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion is characterized by at least one sulfonate group (-S(=O)(=O)[O-])
    attached to an alkane chain (carbon chain).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sulfonate group (-S(=O)(=O)[O-])
    sulfonate_pattern = Chem.MolFromSmarts("[S](=[O])(=[O])[O-]")
    sulfonate_matches = mol.GetSubstructMatches(sulfonate_pattern)
    if not sulfonate_matches:
        return False, "No sulfonate group found"
    
    
    for sulfonate_match in sulfonate_matches:
        sulfur_atom_idx = sulfonate_match[0]
        sulfur_atom = mol.GetAtomWithIdx(sulfur_atom_idx)

    # Check for carbon directly attached to sulfonate sulfur
        carbon_found = False
        for neighbor in sulfur_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                # Check if carbon is aliphatic (sp3 hybridized)
                if neighbor.GetHybridization() == Chem.HybridizationType.SP3:
                    carbon_found = True
                    break #exit the loop when an aliphatic carbon is found

        if not carbon_found:
            return False, "No aliphatic carbon directly attached to sulfonate sulfur"

    return True, "Molecule contains at least one sulfonate group attached to an alkane chain"