"""
Classifies: CHEBI:28892 ganglioside
"""
"""
Classifies: ganglioside
"""

from rdkit import Chem

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    A ganglioside is composed of a glycosphingolipid (ceramide and oligosaccharide) with one or more sialic acids linked on the sugar chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ceramide substructure
    # Ceramide: sphingosine backbone linked via amide bond to fatty acid chain
    # Define a simplified ceramide pattern
    ceramide_smarts = '[NX3][CX3](=O)[C@@H]([C@H](O)[C@H](O)CO)[C@H](O)C=C'
    ceramide_pattern = Chem.MolFromSmarts(ceramide_smarts)
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide substructure found"
    
    # Check for sialic acid substructure
    # Sialic acid: N-acetylneuraminic acid (Neu5Ac) backbone
    sialic_acid_smarts = 'O=C([O-])[C@@H]1[C@@H](O)[C@H](O)[C@H](O[C@H]1CO)NC(C)=O'
    sialic_acid_pattern = Chem.MolFromSmarts(sialic_acid_smarts)
    sialic_acid_matches = mol.GetSubstructMatches(sialic_acid_pattern)
    if len(sialic_acid_matches) == 0:
        return False, "No sialic acid substructure found"

    # Check for glycosidic linkage between ceramide and sugar chain
    glycosidic_linkage_smarts = '[C@H](O[C@H])CO'
    glycosidic_linkage_pattern = Chem.MolFromSmarts(glycosidic_linkage_smarts)
    if not mol.HasSubstructMatch(glycosidic_linkage_pattern):
        return False, "No glycosidic linkage found between ceramide and sugar chain"

    return True, "Contains ceramide backbone, glycosidic-linked sugar chain, and sialic acid residue"