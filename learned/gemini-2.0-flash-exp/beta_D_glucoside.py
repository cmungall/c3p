"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a D-glucopyranose ring
    # The carbon atoms marked with chiral configuration @ and @@ ensure we identify D-glucose.
    # The anomeric carbon is specifically marked as C1 to ensure that we locate the correct stereocenter.
    glucose_pattern = Chem.MolFromSmarts("[C@@H]1([CH2]O)[C@@H](O)[C@H](O)[C@H](O)[C@@H](O)O1")
    matches = mol.GetSubstructMatches(glucose_pattern)

    if not matches:
        return False, "No D-glucose ring found."
    
    # Iterate through all matching glucose rings.
    for match in matches:
        # Get anomeric carbon (C1) index
        anomeric_carbon_index = match[0]
    
        anomeric_atom = mol.GetAtomWithIdx(anomeric_carbon_index)

        #Check for the glycosidic bond connected to C1
        glycosidic_bond = None
        for neighbor in anomeric_atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(anomeric_carbon_index, neighbor.GetIdx())
            if bond.GetBeginAtomIdx() == anomeric_carbon_index and bond.GetEndAtomIdx() != match[1] and bond.GetEndAtomIdx() != match[5] and bond.GetEndAtomIdx() != match[2] and bond.GetEndAtomIdx() != match[3] and bond.GetEndAtomIdx() != match[4]:
                glycosidic_bond = bond
                break
            elif bond.GetEndAtomIdx() == anomeric_carbon_index and bond.GetBeginAtomIdx() != match[1] and bond.GetBeginAtomIdx() != match[5] and bond.GetBeginAtomIdx() != match[2] and bond.GetBeginAtomIdx() != match[3] and bond.GetBeginAtomIdx() != match[4]:
                glycosidic_bond = bond
                break
        
        if glycosidic_bond is None:
              continue

        # Check the chirality of the glycosidic bond at C1. The beta configuration of the anomeric carbon means that
        # the glycosidic bond should point upwards relative to the plane of the glucose ring
        if glycosidic_bond.GetBeginAtomIdx() == anomeric_carbon_index:
           
            if mol.GetAtomWithIdx(glycosidic_bond.GetEndAtomIdx()).GetAtomicNum() == 1:
              continue
            if glycosidic_bond.GetBondDir() == Chem.BondDir.ENDUPRIGHT or glycosidic_bond.GetBondDir() == Chem.BondDir.EITHERDOUBLE: 
                
              return True, "beta-D-glucoside detected (glycosidic bond is ENDUPRIGHT or EITHERDOUBLE)."

        elif glycosidic_bond.GetEndAtomIdx() == anomeric_carbon_index:
            if mol.GetAtomWithIdx(glycosidic_bond.GetBeginAtomIdx()).GetAtomicNum() == 1:
              continue
            if glycosidic_bond.GetBondDir() == Chem.BondDir.ENDDOWNRIGHT or glycosidic_bond.GetBondDir() == Chem.BondDir.EITHERDOUBLE:
              return True, "beta-D-glucoside detected (glycosidic bond is ENDDOWNRIGHT or EITHERDOUBLE)."


    return False, "Not a beta-D-glucoside (no glycosidic bond with correct configuration)"