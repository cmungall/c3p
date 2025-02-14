"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: CHEBI:57575 fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    A fatty acyl-CoA is a fatty acid chain attached to Coenzyme A via a thioester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for the core of the CoA moiety (Adenosine diphosphate and pantetheine)
    # This pattern identifies adenosine diphosphate with a pantetheine attached to it
    coa_pattern = Chem.MolFromSmarts("COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12CC(=O)NCCS")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"
    
    # SMARTS pattern for the thioester linkage
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) != 1:
        return False, f"Found {len(thioester_matches)} thioester group(s), need exactly 1"

    # SMARTS pattern for the fatty acid chain (at least 4 carbons)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)

    # Check if there's any chain of at least 4 carbons attached to the carbonyl group in the thioester
    if not fatty_acid_matches:
      return False, "No fatty acid chain detected"

    # Verify the chain is connected to the thioester
    for match in thioester_matches:
        carbonyl_atom_index = match[0]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_atom_index)

        # check for neighboring carbon chains
        attached_carbon_found = False
        for neighbor in carbonyl_atom.GetNeighbors():
          if neighbor.GetAtomicNum() == 6:
                
                # check the length of the chain
                # traverse the carbons until the chain ends or a branch is found
                visited_atoms = set()
                current_atoms = [neighbor]

                chain_length = 0
                while current_atoms:
                  next_atoms = []
                  for current_atom in current_atoms:
                    visited_atoms.add(current_atom.GetIdx())
                    chain_length += 1
                    #Check neighbors
                    for next_atom in current_atom.GetNeighbors():
                      if next_atom.GetAtomicNum() == 6 and next_atom.GetIdx() not in visited_atoms:
                        next_atoms.append(next_atom)
                  current_atoms = next_atoms
                if chain_length >=4:
                   attached_carbon_found = True
                   break
        if not attached_carbon_found:
          return False, "Fatty acid chain not correctly connected to the thioester group"

    # Check molecular weight - fatty acyl-CoAs typically >700 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, "Molecular weight too low for fatty acyl-CoA"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    if c_count < 10:
         return False, "Too few carbons for fatty acyl-CoA"


    return True, "Contains CoA moiety with a fatty acid chain attached via a thioester bond"