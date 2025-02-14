"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
"""
Classifies: CHEBI:17400 3-oxo-Delta(4) steroid
A 3-oxo steroid conjugated to a C=C double bond at the alpha,beta position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for 3-oxo pattern
    oxo_pattern = Chem.MolFromSmarts("[CX3](=O)")
    oxo_match = mol.GetSubstructMatches(oxo_pattern)
    if not oxo_match:
        return False, "No 3-oxo group found"
    
    # Look for Delta(4) pattern (C=C at positions 4 and 5 in a ring)
    delta4_pattern = Chem.MolFromSmarts("[CR2]=C[CR2]")
    delta4_match = mol.GetSubstructMatches(delta4_pattern)
    if not delta4_match:
        return False, "No Delta(4) alkene found"
    
    # Check if the Delta(4) alkene and 3-oxo group are part of the same fused ring system
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) >= 6:  # Only consider rings with at least 6 atoms
            ring_atoms = set(ring)
            if oxo_match[0][0] in ring_atoms and any(delta4_match[0][i] in ring_atoms for i in range(2)):
                break
    else:
        return False, "Delta(4) alkene and 3-oxo group not part of the same fused ring system"
    
    # Check for correct stereochemistry
    # ... (implement stereochemistry checks)
    
    # Check molecular weight (optional)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 500:
        return False, "Molecular weight outside typical range for steroids"
    
    return True, "Contains a 3-oxo group and a Delta(4) alkene in the same fused ring system"