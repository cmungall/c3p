"""
Classifies: CHEBI:35436 D-glucoside
"""
"""
Classifies: D-glucoside
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside is any glucoside in which the glycoside group is derived from D-glucose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-glucoside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define D-glucose (both alpha and beta anomers)
    alpha_D_glucose_smiles = "OC[C@@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O"
    beta_D_glucose_smiles = "OC[C@@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@H]1O"

    alpha_D_glucose = Chem.MolFromSmiles(alpha_D_glucose_smiles)
    beta_D_glucose = Chem.MolFromSmiles(beta_D_glucose_smiles)

    if alpha_D_glucose is None or beta_D_glucose is None:
        return False, "Failed to generate D-glucose templates"

    # Use substructure search to find D-glucose moiety
    matches_alpha = mol.GetSubstructMatches(alpha_D_glucose)
    matches_beta = mol.GetSubstructMatches(beta_D_glucose)

    if not matches_alpha and not matches_beta:
        return False, "No D-glucose moiety found"

    # Check for glycosidic linkage at the anomeric carbon (C1)
    # The anomeric carbon should be connected to an oxygen which is connected to another group
    glycosidic_bond_found = False
    for match in matches_alpha + matches_beta:
        glucose_atoms = set(match)
        anomeric_carbon_idx = match[1]  # C1 of glucose ring in our templates
        anomeric_carbon = mol.GetAtomWithIdx(anomeric_carbon_idx)

        # Find bonds from the anomeric carbon
        bonds = anomeric_carbon.GetBonds()
        for bond in bonds:
            neighbor = bond.GetOtherAtom(anomeric_carbon)
            if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() not in glucose_atoms:
                # Oxygen connected to anomeric carbon and outside glucose ring
                # Check if oxygen is connected to another non-hydrogen atom outside glucose
                oxygen = neighbor
                oxygen_bonds = oxygen.GetBonds()
                for o_bond in oxygen_bonds:
                    o_neighbor = o_bond.GetOtherAtom(oxygen)
                    if o_neighbor.GetIdx() != anomeric_carbon_idx and o_neighbor.GetAtomicNum() != 1:
                        glycosidic_bond_found = True
                        break
                if glycosidic_bond_found:
                    break
        if glycosidic_bond_found:
            break

    if not glycosidic_bond_found:
        return False, "No glycosidic linkage found at anomeric carbon"

    return True, "Contains D-glucose moiety connected via glycosidic bond"