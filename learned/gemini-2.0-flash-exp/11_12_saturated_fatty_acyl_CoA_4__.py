"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 11,12-saturated fatty acyl-CoA(4-) based on its SMILES string.
    This means checking the fatty acyl chain for saturation at position 11-12.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 11,12-saturated fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # CoA substructure check:
    coa_pattern = Chem.MolFromSmarts('CC(C)(COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12)C(O)C(=O)NCCC(=O)NCCS')
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA substructure"
    
    # Find thioester carbon and number the fatty acid chain carbons:
    thioester_carbon = mol.GetSubstructMatch(Chem.MolFromSmarts('SC(=O)'))
    if not thioester_carbon:
      return False, "Could not identify thioester carbon"
    thioester_carbon = thioester_carbon[1] # second atom in the SMARTS is the C in C=O
    
    fatty_acid_carbon_chain = []
    current_atom = mol.GetAtomWithIdx(thioester_carbon)
    
    # follow the chain
    for _ in range(30): # maximum 30 carbons, avoids infinite loops
        connected_carbons = []
        for neighbor in current_atom.GetNeighbors():
          if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in fatty_acid_carbon_chain:
              connected_carbons.append(neighbor)

        if not connected_carbons:
           break # end of fatty acid chain
        if len(connected_carbons) == 1:
          current_atom = connected_carbons[0]
          fatty_acid_carbon_chain.append(current_atom.GetIdx())
        elif len(connected_carbons) > 1:
          #handle branching, take the non hydroxylated
          for carbon in connected_carbons:
            if not carbon.HasSubstructMatch(Chem.MolFromSmarts('[CX4][OX2H]')):
              current_atom = carbon
              fatty_acid_carbon_chain.append(current_atom.GetIdx())
              break
          else:
             break # branched, do not know what to do.
        
    
    if len(fatty_acid_carbon_chain) < 12: # check chain length.
        return False, "Fatty acid chain too short to have a 11-12 bond"

    # check the bond at 11-12 position
    atom11 = mol.GetAtomWithIdx(fatty_acid_carbon_chain[10]) # 11th carbon of the acyl chain
    atom12 = mol.GetAtomWithIdx(fatty_acid_carbon_chain[11]) # 12th carbon of the acyl chain
    bond = mol.GetBondBetweenAtoms(atom11.GetIdx(), atom12.GetIdx())

    if bond is None:
        return False, "Could not identify bond between carbons 11 and 12"
    if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
        return False, "Bond between carbons 11 and 12 is not a single bond"
    
    return True, "Fatty acyl chain has a saturated bond between carbons 11 and 12 and is a CoA derivative"