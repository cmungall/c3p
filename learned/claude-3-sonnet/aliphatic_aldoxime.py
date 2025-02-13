"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
"""
Classifies: aliphatic aldoxime
Definition: Any aldoxime derived from an aliphatic aldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an aliphatic aldoxime, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First check for aldoxime group pattern [H]C=NOH
    # More specific SMARTS pattern for aldoxime
    aldoxime_pattern = Chem.MolFromSmarts("[CH1](=N[OH1])-[!$(C=O);!$(C=N);!$(C=C);!$([N+]=[O-])]")
    
    if not mol.HasSubstructMatch(aldoxime_pattern):
        return False, "No proper aldoxime group (HC=NOH) found"

    # Get matches for the aldoxime group
    matches = mol.GetSubstructMatches(aldoxime_pattern)
    
    # Check each potential aldoxime group
    for match in matches:
        c_atom = mol.GetAtomWithIdx(match[0])
        n_atom = mol.GetAtomWithIdx(match[1])
        o_atom = mol.GetAtomWithIdx(match[2])
        
        # Verify hybridization
        if c_atom.GetHybridization() != Chem.HybridizationType.SP2:
            return False, "Carbon in C=N must be sp2 hybridized"
        
        # Check that carbon has exactly one hydrogen
        if c_atom.GetTotalNumHs() != 1:
            return False, "Carbon must have exactly one hydrogen (aldoxime)"
        
        # Check formal charges
        if c_atom.GetFormalCharge() != 0 or n_atom.GetFormalCharge() != 0:
            return False, "Unexpected formal charges on C or N atoms"
            
        # Check for aci-nitro compounds
        aci_nitro = Chem.MolFromSmarts("[O-][N+](=C)O")
        if mol.HasSubstructMatch(aci_nitro):
            return False, "Compound is an aci-nitro compound, not an aldoxime"

    # Check for aromatic or conjugated systems
    conjugated_patterns = [
        Chem.MolFromSmarts("C=C-C=N"), # conjugated alkene
        Chem.MolFromSmarts("C=C-N=C"), # conjugated imine
        Chem.MolFromSmarts("O=C-C=N"), # conjugated carbonyl
        Chem.MolFromSmarts("n1ccccc1"), # pyridine
        Chem.MolFromSmarts("c1ccccc1")  # benzene
    ]
    
    for pattern in conjugated_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return False, "Contains conjugated or aromatic system"

    # Check that carbon attached to C=N is sp3 (aliphatic)
    for match in matches:
        neighbors = mol.GetAtomWithIdx(match[0]).GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetIdx() not in [match[1]]:  # exclude the N atom
                if neighbor.GetHybridization() != Chem.HybridizationType.SP3:
                    return False, "Carbon chain must be aliphatic"

    return True, "Contains aliphatic aldoxime group (HC=NOH)"