"""
Classifies: CHEBI:28863 flavanones
"""
"""
Classifies: CHEBI:27555 flavanone
"""

from rdkit import Chem
from rdkit.Chem import rdFMCS

def is_flavanone(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    Flavanones are defined as members of the class of flavans with a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the core flavanone scaffold
    flavanone_core = Chem.MolFromSmarts("[C@@]12(C=CC=C1)OC(=O)C(C2)c1ccccc1")
    
    # Check if the molecule contains the flavanone core
    if mol.HasSubstructMatch(flavanone_core):
        return True, "Contains the 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one flavanone skeleton"
    else:
        # Check for common substituted derivatives
        derivatives = [
            Chem.MolFromSmarts("[C@]12(C=CC=C1)OC(=O)C(C2)c1ccccc1"), # 2-aryl-2H-chromen-4-one (flavanone with additional chiral center)
            Chem.MolFromSmarts("[C@@]12(C=CC=C1)OC(=O)C(C2)c1ccccc1O"), # 7-hydroxyflavanone
            Chem.MolFromSmarts("[C@@]12(C=CC=C1)OC(=O)C(C2)c1ccc(O)cc1"), # 4'-hydroxyflavanone
            # Add more common derivative patterns here
        ]
        
        for derivative in derivatives:
            if mol.HasSubstructMatch(derivative):
                return True, "Contains a substituted flavanone derivative"
        
        return False, "Does not contain the flavanone skeleton or common derivatives"