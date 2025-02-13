"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: CHEBI:36737 quinone

A quinone is a compound having a fully conjugated cyclic dione structure, such as that of 
benzoquinones, derived from aromatic compounds by conversion of an even number of -CH= groups 
into -C(=O)- groups with any necessary rearrangement of double bonds (polycyclic and 
heterocyclic analogues are included).
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdqueries

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for conjugated cyclic dione pattern
    dione_pattern = rdqueries.GetQueryForSubstructureMatch(Chem.MolFromSmarts("[C&R1]=,:[C&R1]=,:[C&R1]=,:[C&R1]=,:[C&R1]=,:[C&R1]=,:[C&R1]=,:[C&R1]=,:1]=[O&R1]:[a&r]:[O&R1]:[a&r]1"))
    if not mol.HasSubstructMatch(dione_pattern):
        return False, "No conjugated cyclic dione structure found"
    
    # Count aromatic rings
    num_aromatic_rings = AllChem.CalcNumAromaticRings(mol)
    if num_aromatic_rings == 0:
        return False, "No aromatic rings present"
    
    # Check if carbonyl groups are part of an aromatic ring
    carboxyl_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0]
    for atom in carboxyl_atoms:
        if not atom.IsInRingOfSize(6):
            return False, "Carbonyl groups not part of an aromatic ring"
    
    return True, "Contains a fully conjugated cyclic dione structure derived from an aromatic compound"