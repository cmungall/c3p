"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
"""
Classifies: CHEBI:27732 beta-D-glucosiduronic acid
A glucosiduronic acid resulting from the formal condensation of any substance 
with beta-D-glucuronic acid to form a glycosidic bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for beta-D-glucuronic acid substructure
    glucuronic_pattern = Chem.MolFromSmarts("[C@@H]1([C@@H](O)[C@H](O)[C@@H](O)[C@H](O1)O)[C@H](O)C(O)=O")
    if not mol.HasSubstructMatch(glucuronic_pattern):
        return False, "No beta-D-glucuronic acid substructure found"
    
    # Look for glycosidic bond (O-C)
    glycosidic_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if not glycosidic_matches:
        return False, "No glycosidic bond found"
    
    # Check for carbon skeleton attached to glycosidic oxygen
    for match in glycosidic_matches:
        o_atom = mol.GetAtomWithIdx(match[0])
        c_atom = mol.GetAtomWithIdx(match[1])
        if c_atom.GetAtomicNum() == 6 and c_atom.GetTotalNumHs() <= 1:
            return True, "Contains beta-D-glucuronic acid substructure linked via glycosidic bond"
    
    return False, "No carbon skeleton attached to glycosidic oxygen"