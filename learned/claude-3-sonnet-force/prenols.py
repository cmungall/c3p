"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: CHEBI:38732 prenols
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    Prenols are alcohols with the general formula H-[CH2C(Me)=CHCH2]nOH,
    where the carbon skeleton is composed of one or more isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OX1H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"
    
    # Look for isoprene unit pattern (CH2=C(CH3)CH=CH2)
    isoprene_pattern = Chem.MolFromSmarts("[CH2]=[C]([CH3])[CH]=[CH2]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) == 0:
        return False, "No isoprene units found"
    
    # Check that isoprene units are connected in a linear chain
    isoprene_atoms = set([atom.GetIdx() for match in isoprene_matches for atom in match])
    isoprene_bonds = [bond for bond in mol.GetBonds() if bond.GetBeginAtomIdx() in isoprene_atoms and bond.GetEndAtomIdx() in isoprene_atoms]
    if len(isoprene_bonds) != len(isoprene_matches) * 2:
        return False, "Isoprene units not connected in a linear chain"
    
    # Check for hydroxyl group at one end of the isoprene chain
    start_atom = mol.GetAtomWithIdx(isoprene_matches[0][0])
    end_atom = mol.GetAtomWithIdx(isoprene_matches[-1][-1])
    has_start_hydroxyl = start_atom.GetTotalNumHs() == 1 and sum(bond.GetBondTypeAsDouble() for bond in mol.GetBondEdges(start_atom.GetIdx())) == 2
    has_end_hydroxyl = end_atom.GetTotalNumHs() == 1 and sum(bond.GetBondTypeAsDouble() for bond in mol.GetBondEdges(end_atom.GetIdx())) == 2
    if not (has_start_hydroxyl or has_end_hydroxyl):
        return False, "Hydroxyl group not at end of isoprene chain"
    
    return True, "Molecule contains a linear chain of isoprene units with a terminal hydroxyl group"