"""
Classifies: CHEBI:39437 tocol
"""
"""
Classifies: tocol compounds
Definition: A chromanol with a chroman-6-ol skeleton substituted at position 2 
by a saturated or triply-unsaturated hydrocarbon chain of three isoprenoid units
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_tocol, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for basic chroman core (benzene fused to pyran)
    chroman_pattern = Chem.MolFromSmarts("O1CCCc2ccccc12")
    if not mol.HasSubstructMatch(chroman_pattern):
        return False, "No chroman core found"
    
    # Check for hydroxyl group (or modified hydroxyl) at position 6
    # Note: Position 6 corresponds to para position relative to oxygen
    hydroxy_pattern = Chem.MolFromSmarts("O1CCCc2cc([OH,O])ccc12")
    hydroxy_ester_pattern = Chem.MolFromSmarts("O1CCCc2cc(O[C,P])ccc12")
    if not (mol.HasSubstructMatch(hydroxy_pattern) or mol.HasSubstructMatch(hydroxy_ester_pattern)):
        return False, "No hydroxyl group (or derivative) at position 6"
    
    # Check for substitution at position 2 (carbon next to oxygen in pyran ring)
    # Looking for branched carbon chain
    c2_substitution = Chem.MolFromSmarts("O1[C](CC)CCc2ccccc12")
    if not mol.HasSubstructMatch(c2_substitution):
        return False, "No substitution at position 2"
    
    # Count carbons in the longest chain from position 2
    # This checks for the presence of the isoprenoid chain
    side_chain_pattern = Chem.MolFromSmarts("[CH3]C(C)CCC[CH](C)CCC[CH](C)CCC")
    side_chain_unsat = Chem.MolFromSmarts("[CH3]C(C)=CCC[CH](C)=CCC[CH](C)=CCC")
    
    if not (mol.HasSubstructMatch(side_chain_pattern) or mol.HasSubstructMatch(side_chain_unsat)):
        return False, "Missing required isoprenoid chain"
    
    # Check total number of carbons to ensure proper chain length
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 20:  # Minimum carbons for a tocol
        return False, "Insufficient carbons for tocol structure"
        
    # Additional check for common modifications that don't invalidate tocol status
    valid_modifications = [
        ("acetate", "CC(=O)O"),
        ("succinate", "O=C(O)CCC(=O)O"),
        ("carboxyl", "C(=O)O"),
    ]
    
    # Look for common modifications
    modifications = []
    for mod_name, mod_smarts in valid_modifications:
        pattern = Chem.MolFromSmarts(mod_smarts)
        if mol.HasSubstructMatch(pattern):
            modifications.append(mod_name)
    
    base_message = "Contains chroman-6-ol core with appropriate isoprenoid chain"
    if modifications:
        return True, f"{base_message} (modified with: {', '.join(modifications)})"
    
    return True, base_message