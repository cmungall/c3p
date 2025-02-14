"""
Classifies: CHEBI:26377 pterocarpans
"""
"""
Classifies: CHEBI:26395 pterocarpans
"""
from rdkit import Chem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    A pterocarpan is characterized by a 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pterocarpan, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for pterocarpan core
    # The pattern represents the fused tetracyclic ring system
    pterocarpan_smarts = """
    [#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1
    -2-[#8]-[#6]3-[#6]-[#6]-[#8]-[#6]-4-[#6]:[#6]:[#6]:[#6]:[#6]:[#6]:4-[#6]-3-2
    """

    # Remove whitespace and line breaks
    pterocarpan_smarts = pterocarpan_smarts.replace('\n', '').replace(' ', '')
    
    pterocarpan_core = Chem.MolFromSmarts(pterocarpan_smarts)
    if pterocarpan_core is None:
        return False, "Could not define pterocarpan core pattern"

    # Check for substructure match
    if mol.HasSubstructMatch(pterocarpan_core):
        return True, "Molecule contains the pterocarpan core structure"

    return False, "Molecule does not contain the pterocarpan core structure"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:26395',
        'name': 'pterocarpans',
        'definition': 'Members of the class of benzofurochromene with a 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton and its substituted derivatives. They generally bear structural resemblance to isoflavanoids that possess antibiotic activity and are produced by plant tissues in response to infection. They are the 3,4-dihydroderivatives of coumestans.',
        'parents': []
    }
}