"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA is an unsaturated fatty acyl-CoA with a double bond between positions 2 and 3
    relative to the carbonyl of the thioester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the CoA substructure (simplified, focusing on the key parts, Adenosine removed to avoid issues)
    coa_pattern = Chem.MolFromSmarts("CC(C)(COP(=O)([OX1])OP(=O)([OX1])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)([OX1])[OX1])n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCS")
    
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA substructure not found"

    # Locate the thioester group (-C(=O)-S-) to begin fatty acid chain analysis
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
      return False, "No thioester group found"

    # Check for the alpha-beta unsaturation (C=C-C=O)
    alpha_beta_pattern = Chem.MolFromSmarts("[C]=[C]-[C](=[O])")
    alpha_beta_matches = mol.GetSubstructMatches(alpha_beta_pattern)
    if not alpha_beta_matches:
      return False, "No double bond between position 2 and 3 (alpha/beta) of thioester"


    #check for at least 4 carbons after the double bond, ensure a plausible chain length
    acyl_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]=[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    acyl_chain_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if not acyl_chain_matches:
        return False, "Fatty acid chain seems too short or missing"


    # Additional check: ensure the double bond is directly attached to the C of the carbonyl
    # First get the C of the carbonyl for the identified thioester match:
    for match in thioester_matches:
        carbonyl_c_idx = match[0]
        # now get the atoms attached to it:
        neighbors = mol.GetAtomWithIdx(carbonyl_c_idx).GetNeighbors()
        
        #check if there is at least one carbon, and if it has a double bond
        found_double_bond = False
        for n in neighbors:
          if n.GetAtomicNum() == 6 and  mol.GetBondBetweenAtoms(carbonyl_c_idx, n.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
            found_double_bond = True
            break
        if not found_double_bond:
            return False, "Double bond is not directly attached to carbonyl group"




    return True, "Molecule contains CoA, a thioester, and a double bond between the 2nd and 3rd carbon atoms of the acyl chain."