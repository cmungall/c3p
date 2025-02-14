"""
Classifies: CHEBI:64985 bioconjugate
"""
"""
Classifies: CHEBI:36357 bioconjugate
A molecular entity consisting of at least 2 biological molecules covalently linked together.
"""
from rdkit import Chem
from rdkit.Chem import FragmentMatcher
from rdkit.Chem import rdFMCS

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate is a molecular entity consisting of at least two biological molecules covalently linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get fragments
    fragments = list(Chem.rdmolops.GetMolFrags(mol, aromatics=Chem.rdchem.AromaticityType.SIMPLE))
    
    # Count number of "biological" fragments
    biological_fragments = []
    for frag in fragments:
        frag_mol = Chem.rdmolops.MolFragmentToSmiles(mol, frag, kekuleSmiles=True)
        frag_mol = Chem.MolFromSmiles(frag_mol)
        
        # Check if fragment contains common biological substructures
        biological_patterns = ['[NX3]', '[NX4+]', 'C(=O)[O,N]', 'n', 'O=C-O-P', 'S']
        is_biological = any(frag_mol.HasSubstructMatch(Chem.MolFromSmarts(patt)) for patt in biological_patterns)
        
        if is_biological:
            biological_fragments.append(frag_mol)
    
    # Check if there are at least 2 biological fragments
    if len(biological_fragments) < 2:
        return False, "Less than 2 biological fragments found"
    
    # Check if biological fragments are linked
    fmcs = rdFMCS.FindMCS(biological_fragments, matchValences=True, completeRingsOnly=True)
    if fmcs.numBonds == 0:
        return False, "Biological fragments are not covalently linked"
    
    return True, "Contains at least 2 biological molecules covalently linked together"