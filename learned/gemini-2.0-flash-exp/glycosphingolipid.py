"""
Classifies: CHEBI:24402 glycosphingolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    A glycosphingolipid is a glycolipid that is a carbohydrate-containing derivative of a sphingoid or ceramide.
    The carbohydrate residue is attached by a glycosidic linkage to O-1 of the sphingoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosphingolipid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for sphingosine/ceramide core
    # This pattern checks for a long carbon chain, an amide bond, a hydroxyl group, and a linkage to an oxygen-bearing carbon
    # The hydroxyl group on C2 is not strictly necessary in a sphingosine, but many glycosphingolipids are like that.
    # Also added an extra C=C bond to account for sphingosine
    ceramide_pattern = Chem.MolFromSmarts("[CX4][CX4](O)[CX4]([NX3][CX3](=[OX1])[#6])([CX4][O][CX4])([#6])=[#6]") #Added a final double bond, sphingosine
    if not mol.HasSubstructMatch(ceramide_pattern):
        ceramide_pattern_nodb = Chem.MolFromSmarts("[CX4][CX4](O)[CX4]([NX3][CX3](=[OX1])[#6])([CX4][O][CX4])([#6])")
        if not mol.HasSubstructMatch(ceramide_pattern_nodb):
            return False, "No ceramide/sphingosine core found"

    # SMARTS pattern for glycosidic linkage
    # This pattern looks for a carbon attached to an oxygen, where the oxygen is also part of a carbohydrate ring
    glycosidic_pattern = Chem.MolFromSmarts("[OX2][CX4]([OX2])[CX4]1([OX2][CX4]([OX2])[CX4][CX4]([OX2])[CX4]1)") #Modified the pattern to detect a sugar ring
    if not mol.HasSubstructMatch(glycosidic_pattern):
         return False, "No glycosidic linkage found"
    
    # Check for long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be sphingolipid"

    # Check molecular weight - glycosphingolipids typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
         return False, "Molecular weight too low for a glycosphingolipid"


    return True, "Contains a ceramide/sphingosine core with a carbohydrate attached via a glycosidic linkage"