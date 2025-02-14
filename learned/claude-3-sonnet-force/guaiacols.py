"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: CHEBI:51504 guaiacols
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_guaiacol(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol is any phenol carrying an additional methoxy substituent at the ortho-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for phenol ring pattern (aromatic ring with -OH group)
    phenol_pattern = Chem.MolFromSmarts("c1ccccc1O")
    if not mol.HasSubstructMatch(phenol_pattern):
        return False, "Not a phenol"
    
    # Look for methoxy group (-O-C) at ortho position
    methoxy_pattern = Chem.MolFromSmarts("Oc1ccccc1OC")
    if not mol.HasSubstructMatch(methoxy_pattern):
        return False, "No ortho-methoxy group found"
    
    # Check if there is only one methoxy group
    methoxy_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("OC"))
    if len(methoxy_groups) != 1:
        return False, f"Found {len(methoxy_groups)} methoxy groups, should be exactly 1"
    
    # Check if the methoxy group is ortho to the -OH
    methoxy_atom = mol.GetAtoms()[methoxy_groups[0][0]].GetIdx()
    hydroxy_atom = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "O" and atom.GetFormalCharge() == 0][0]
    
    if not mol.GetBondBetweenAtoms(methoxy_atom, hydroxy_atom).IsInRingOfSize(6):
        return False, "Methoxy group not ortho to phenolic -OH"
    
    return True, "Molecule contains a phenol with an ortho-methoxy substituent"