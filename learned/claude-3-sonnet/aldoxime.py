"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: CHEBI:33566 aldoxime
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    Aldoximes are compounds with the general structure R-CH=N-OH, where R is any carbon group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for oxime group pattern (C=N-OH)
    # [CX3H1] ensures carbon has exactly one hydrogen
    # [NX2] ensures nitrogen has double bond
    # [OX2H1] ensures oxygen is connected to one hydrogen
    oxime_pattern = Chem.MolFromSmarts("[CX3H1]=[NX2]-[OX2H1]")
    
    # Alternative pattern for deprotonated oximes
    oxime_anion_pattern = Chem.MolFromSmarts("[CX3H1]=[NX2]-[OX1-]")
    
    matches = mol.GetSubstructMatches(oxime_pattern)
    anion_matches = mol.GetSubstructMatches(oxime_anion_pattern)
    
    total_matches = len(matches) + len(anion_matches)
    
    if total_matches == 0:
        return False, "No aldoxime group (C(H)=N-OH) found"
    
    # Check if the carbon of C=N is connected to exactly one hydrogen
    # This distinguishes aldoximes from ketoximes
    for match in matches:
        carbon_atom = mol.GetAtomWithIdx(match[0])
        if carbon_atom.GetTotalNumHs() != 1:
            return False, "Oxime carbon has incorrect number of hydrogens"
            
    for match in anion_matches:
        carbon_atom = mol.GetAtomWithIdx(match[0])
        if carbon_atom.GetTotalNumHs() != 1:
            return False, "Oxime carbon has incorrect number of hydrogens"
    
    # Verify carbon is not part of C=O group (would indicate an oxime ester)
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=O)-[OX2]")
    if mol.HasSubstructMatch(carbonyl_pattern):
        # Check if the carbonyl carbon is the same as oxime carbon
        carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
        for oxime_match in matches:
            oxime_carbon = oxime_match[0]
            for carbonyl_match in carbonyl_matches:
                if oxime_carbon == carbonyl_match[0]:
                    return False, "Contains oxime ester group instead of aldoxime"
    
    return True, "Contains aldoxime group (R-CH=N-OH)"