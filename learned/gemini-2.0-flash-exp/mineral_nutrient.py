"""
Classifies: CHEBI:167164 mineral nutrient
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    Mineral nutrients are inorganic compounds that provide essential elements for the body.
    They generally contain metals or other charged inorganic compounds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral nutrient, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define common mineral nutrient metals (alkali, alkaline earth, common transition metals)
    mineral_metals = [3, 11, 12, 19, 20, 25, 26, 29, 30, 37, 38, 55, 56] # Li, Na, Mg, K, Ca, Mn, Fe, Cu, Zn, Rb, Sr, Cs, Ba

    #Check if metals are present
    metal_atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() in mineral_metals]
    
    #If metals are absent, check if common counterions (halides, etc) and common minerals are present.
    if not metal_atoms:
        
        halide_pattern = Chem.MolFromSmarts("[F-,Cl-,Br-,I-]")
        sulfate_pattern = Chem.MolFromSmarts("[O-]S([O-])(=O)=O")
        phosphate_pattern = Chem.MolFromSmarts("[O-]P([O-])([O-])=O")
        carbonate_pattern = Chem.MolFromSmarts("[O-]C([O-])=O")
        nitrate_pattern = Chem.MolFromSmarts("[O-][N+]([O-])=O")
        silicate_pattern = Chem.MolFromSmarts("[O-][Si]([O-])([O-])[O-]")
        hydroxide_pattern = Chem.MolFromSmarts("[OH-]")
        oxide_pattern = Chem.MolFromSmarts("[O--]")
        sulfide_pattern = Chem.MolFromSmarts("[S--]")

        if not (mol.HasSubstructMatch(halide_pattern) or mol.HasSubstructMatch(sulfate_pattern) or
                mol.HasSubstructMatch(phosphate_pattern) or mol.HasSubstructMatch(carbonate_pattern) or
                mol.HasSubstructMatch(nitrate_pattern) or mol.HasSubstructMatch(silicate_pattern) or
                mol.HasSubstructMatch(hydroxide_pattern) or mol.HasSubstructMatch(oxide_pattern) or
                mol.HasSubstructMatch(sulfide_pattern) ):

           return False, "No metal, or common mineral ions found"
   
    #Check for 'small' organic fragments (less than 5 non-H atoms)
    organic_pattern = Chem.MolFromSmarts("[#6]") #check for carbons to see if there is an organic part
    if mol.HasSubstructMatch(organic_pattern):
         for frag in Chem.GetMolFrags(mol, asMols=True):
            if frag.HasSubstructMatch(Chem.MolFromSmarts("[#6]")): #If this fragment is organic
               num_non_h_atoms = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() != 1)
               if num_non_h_atoms > 5: #If this organic fragment is too large
                 return False, "Contains large organic fragments, not typical for a mineral nutrient"

    # Check for charge balance of the entire molecule.
    total_charge = Chem.GetFormalCharge(mol)
    if total_charge != 0:
        return False, "Overall charge is not neutral"


    return True, "Contains metal ions and/or typical mineral ions"