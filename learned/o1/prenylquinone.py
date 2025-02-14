"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: prenylquinone
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone is a quinone substituted by a polyprenyl-derived side-chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define quinone core patterns (benzoquinone and naphthoquinone)
    benzoquinone_pattern = Chem.MolFromSmarts('O=C1C=CC(=O)C=C1')  # 1,4-benzoquinone
    naphthoquinone_pattern = Chem.MolFromSmarts('O=C1C=CC2=CC=CC(=O)C2=C1')  # 1,4-naphthoquinone
    
    # Check for the presence of quinone core
    if mol.HasSubstructMatch(benzoquinone_pattern):
        quinone_type = 'benzoquinone'
        quinone_match = mol.GetSubstructMatch(benzoquinone_pattern)
    elif mol.HasSubstructMatch(naphthoquinone_pattern):
        quinone_type = 'naphthoquinone'
        quinone_match = mol.GetSubstructMatch(naphthoquinone_pattern)
    else:
        return False, "No quinone core (benzoquinone or naphthoquinone) found"
    
    # Get the atoms in the quinone core
    quinone_atoms = set(quinone_match)
    
    # Identify side chains attached to quinone core
    # Find attachment points
    attachments = []
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtomIdx()
        end_atom = bond.GetEndAtomIdx()
        if (begin_atom in quinone_atoms and end_atom not in quinone_atoms) or \
           (end_atom in quinone_atoms and begin_atom not in quinone_atoms):
            # This bond connects the quinone core to the rest of the molecule
            attachments.append(bond)
    
    # Check if any of the attachments lead to a polyprenyl chain
    has_polyprenyl_chain = False
    for bond in attachments:
        # Get the side chain starting atom
        side_chain_atom_idx = bond.GetBeginAtomIdx() if bond.GetEndAtomIdx() in quinone_atoms else bond.GetEndAtomIdx()
        
        # Use a depth-limited search to find polyprenyl chain patterns
        side_chain = Chem.PathToSubmol(mol, Chem.FindAtomEnvironmentOfRadiusN(mol, radius=20, atomIdx=side_chain_atom_idx))
        
        # Define isoprene unit pattern (C=C-C-C=C)
        isoprene_unit_pattern = Chem.MolFromSmarts('C(=C)C-C=C')
        
        # Count the number of isoprene units in the side chain
        isoprene_matches = side_chain.GetSubstructMatches(isoprene_unit_pattern)
        if len(isoprene_matches) >= 2:
            has_polyprenyl_chain = True
            break
    
    if not has_polyprenyl_chain:
        return False, "No polyprenyl-derived side chain found attached to the quinone core"
    
    return True, f"Contains {quinone_type} core with polyprenyl-derived side chain"