"""
Classifies: CHEBI:26004 phenylpropanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Core Phenylpropane Structure - more flexible now
    phenylpropane_pattern1 = Chem.MolFromSmarts("c1ccccc1[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    phenylpropane_pattern2 = Chem.MolFromSmarts("c1ccccc1[CX4,CX3]~[CX4,CX3]")


    # 2. Common Modifications (hydroxyl, carbonyl, methoxy, ester etc)
    modification_pattern1 = Chem.MolFromSmarts("[OH]")
    modification_pattern2 = Chem.MolFromSmarts("[CX3]=[OX1]")
    modification_pattern3 = Chem.MolFromSmarts("[CX4]-[OX2]")
    modification_pattern4 = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    modification_pattern5 = Chem.MolFromSmarts("O[C@H]1[C@@H]([C@@H](O)[C@H](O)[C@H](CO)O)O[C@H]1[OX2]") # glycoside


    # 3. Common ring systems:
    benzopyran_pattern = Chem.MolFromSmarts("c1ccccc1-C=2OC=CC2") # more specific
    benzodioxole_pattern = Chem.MolFromSmarts("c1ccccc1-C2-O-C-O-2")
    isobenzofuranone_pattern = Chem.MolFromSmarts("c1ccccc1-C2-C(=O)-O-2")
    coumarin_pattern = Chem.MolFromSmarts("c1ccccc1-C=2C=CC(=O)O2")
    flavone_pattern = Chem.MolFromSmarts("c1ccccc1-c2oc(cc(=O)c2c1)-c3ccccc3")


    has_core = False
    if mol.HasSubstructMatch(phenylpropane_pattern1) or mol.HasSubstructMatch(phenylpropane_pattern2):
        has_core = True


    if has_core:
       if (modification_pattern1 and mol.HasSubstructMatch(modification_pattern1)) or \
           (modification_pattern2 and mol.HasSubstructMatch(modification_pattern2)) or \
            (modification_pattern3 and mol.HasSubstructMatch(modification_pattern3)) or \
            (modification_pattern4 and mol.HasSubstructMatch(modification_pattern4)) or \
            (modification_pattern5 and mol.HasSubstructMatch(modification_pattern5)) or \
            (benzopyran_pattern and mol.HasSubstructMatch(benzopyran_pattern)) or \
            (benzodioxole_pattern and mol.HasSubstructMatch(benzodioxole_pattern)) or \
            (isobenzofuranone_pattern and mol.HasSubstructMatch(isobenzofuranone_pattern)) or \
             (coumarin_pattern and mol.HasSubstructMatch(coumarin_pattern)) or \
             (flavone_pattern and mol.HasSubstructMatch(flavone_pattern)) :
           return True, "Phenylpropane core with typical modification or heterocyclic ring system detected."
       else:
            return True, "Phenylpropane core detected, but no explicit phenylpropanoid modification found, but core suffices."
    
    
    return False, "No phenylpropane core (C6-C2 or C6-C3) structure detected."