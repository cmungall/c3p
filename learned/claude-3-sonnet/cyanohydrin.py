"""
Classifies: CHEBI:23437 cyanohydrin
"""
"""
Classifies: CHEBI:35667 cyanohydrin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cyanohydrin(smiles: str):
    """
    Determines if a molecule is a cyanohydrin based on its SMILES string.
    A cyanohydrin is an alpha-hydroxynitrile resulting from the formal addition 
    of hydrogen cyanide to the C=O bond of an aldehyde or ketone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyanohydrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carbon with both OH and CN groups
    # [C;X4] = sp3 carbon (any substitution)
    # [OX2H1] = hydroxy group (-OH)
    # [CX2]#[NX1] = cyano group (-C≡N)
    # The OH and CN should be on the same carbon but with any other substituents
    cyanohydrin_pattern = Chem.MolFromSmarts('[C;X4;!$(C1[O,N]CC1)]([OX2H1])-[CX2]#[NX1]')
    
    if not mol.HasSubstructMatch(cyanohydrin_pattern):
        return False, "No alpha-hydroxynitrile group found"
    
    # Count number of such patterns - should be exactly one
    matches = mol.GetSubstructMatches(cyanohydrin_pattern)
    if len(matches) > 1:
        return False, "Multiple alpha-hydroxynitrile groups found"
    
    # Verify carbon hybridization of the central carbon
    match = matches[0]
    c_atom = mol.GetAtomWithIdx(match[0])
    if c_atom.GetHybridization() != Chem.HybridizationType.SP3:
        return False, "Carbon with OH and CN must be sp3 hybridized"
        
    # Check that molecule is not too small (should be larger than HCN)
    if mol.GetNumAtoms() < 4:
        return False, "Molecule too small to be a cyanohydrin"

    # Additional check to exclude cyclic structures where OH and CN are part of a small ring
    ring_pattern = Chem.MolFromSmarts('[C;R]1[O,N][C;R]1')
    if mol.HasSubstructMatch(ring_pattern):
        return False, "OH and CN groups are part of a ring system"
        
    return True, "Contains alpha-hydroxynitrile group derived from aldehyde/ketone"