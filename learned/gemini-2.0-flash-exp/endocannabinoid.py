"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is likely an endocannabinoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely an endocannabinoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for long carbon chain (>= 10C) and unsaturation. Use # rotatable bonds as a proxy.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8: #adjust this value
        return False, "Too few rotatable bonds (likely short chain)"
    
    # Pattern to find long carbon chains with double bonds
    long_chain_pattern = Chem.MolFromSmarts("C~C=C~C") #C-C=C-C or longer
    if not mol.HasSubstructMatch(long_chain_pattern):
         long_chain_pattern = Chem.MolFromSmarts("C~C~C~C~C~C~C~C~C~C") #also checking for very long chains without unsat.
         if not mol.HasSubstructMatch(long_chain_pattern):
              return False, "No long carbon chain with or without unsaturation found"
    

    # 2. Check for polar head group (ethanolamide or glycerol) and linkage (ester or amide).
    ethanolamide_pattern = Chem.MolFromSmarts("NCCO")
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])") # changed from NX3 to allow N w/ two bonds to be recognized

    has_ethanolamide = mol.HasSubstructMatch(ethanolamide_pattern)
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_amide = mol.HasSubstructMatch(amide_pattern)

    #Check connection to long carbon chain
    if has_amide:
       amide_matches = mol.GetSubstructMatches(amide_pattern)
       for match in amide_matches:
          carbon_atoms_linked_to_amide_N = [mol.GetAtomWithIdx(neighbor).GetAtomicNum() for neighbor in mol.GetAtomWithIdx(match[0]).GetNeighbors() if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6]
          if len(carbon_atoms_linked_to_amide_N) > 0 :
            #found a carbon directly attached
             break
       else:
           return False, "No long chain attached to amide nitrogen" # no carbon directly connected to amide nitrogen
    if has_ester:
       ester_matches = mol.GetSubstructMatches(ester_pattern)
       for match in ester_matches:
            carbon_atoms_linked_to_ester_O = [mol.GetAtomWithIdx(neighbor).GetAtomicNum() for neighbor in mol.GetAtomWithIdx(match[0]).GetNeighbors() if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6]
            if len(carbon_atoms_linked_to_ester_O) > 0:
                break
       else:
            return False, "No long chain attached to ester oxygen"


    if not (has_ethanolamide or has_glycerol) or not (has_ester or has_amide):
        return False, "No relevant polar head group or linkage found"


    # 3. If both conditions met, its likely an endocannabinoid
    return True, "Likely an endocannabinoid based on long chain and polar head group"