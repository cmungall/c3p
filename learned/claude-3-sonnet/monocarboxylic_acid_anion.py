"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
"""
Classifies: CHEBI:36342 monocarboxylic acid anion

A carboxylic acid anion formed when the carboxy group of a monocarboxylic acid is deprotonated.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylate anion group (-C(=O)[O-])
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2-,OX1-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(carboxylate_matches) != 1:
        return False, f"Found {len(carboxylate_matches)} carboxylate groups, need exactly 1"
    
    # Check for the presence of a hydrocarbon chain or ring system
    hydrocarbon_pattern = Chem.MolFromSmarts("[C;!$(C=O)]")
    hydrocarbon_matches = mol.GetSubstructMatches(hydrocarbon_pattern)
    if len(hydrocarbon_matches) < 3:
        return False, "Not enough hydrocarbon atoms to be a monocarboxylic acid anion"
    
    # Check for the presence of hydroxy groups (-OH)
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    # Check for the presence of carbon-carbon double bonds (C=C)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    
    # Check molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 60 or mol_wt > 500:
        return False, "Molecular weight outside typical range for monocarboxylic acid anions"
    
    hba = rdMolDescriptors.CalcNumLipinskiHBA(mol)
    hbd = rdMolDescriptors.CalcNumLipinskiHBD(mol)
    if hba < 2 or hba > 10 or hbd < 1 or hbd > 6:
        return False, "Hydrogen bond donor/acceptor counts outside typical range"
    
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    if tpsa < 20 or tpsa > 150:
        return False, "Topological polar surface area outside typical range"
    
    # If all checks pass, classify as a monocarboxylic acid anion
    reason = "Contains a single carboxylate anion group, a hydrocarbon chain/ring, and molecular properties consistent with a monocarboxylic acid anion."
    if len(hydroxy_matches) > 0:
        reason += " Also contains hydroxy groups."
    if len(double_bond_matches) > 0:
        reason += " Also contains carbon-carbon double bonds."
    return True, reason