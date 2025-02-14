"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    A trans-2-enoyl-CoA is a coenzyme A molecule with a fatty acid chain
    attached via a thioester bond with a trans double bond at the 2-position
    relative to the carbonyl.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trans-2-enoyl-CoA, False otherwise.
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the CoA SMARTS pattern. This is the minimal core of the CoA moiety
    # we are checking to avoid issues with tautomers.
    coa_pattern = Chem.MolFromSmarts('SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12')
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"
    
    # Define the thioester group
    thioester_pattern = Chem.MolFromSmarts("S[CX3](=O)")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) != 1:
         return False, f"Found {len(thioester_matches)} thioester groups, need exactly 1"


    # Define the trans double bond pattern
    trans_double_bond_pattern = Chem.MolFromSmarts("[CX4]=[CX4]")
    
    
    # Find the thioester carbon and check the adjacent atoms
    for match in thioester_matches:
        thioester_carbon_index = match[1]
        thioester_carbon = mol.GetAtomWithIdx(thioester_carbon_index)

        # Find the two atoms connected to the thioester carbon
        neighboring_atoms = [atom.GetIdx() for atom in thioester_carbon.GetNeighbors()]
        if len(neighboring_atoms) != 3:
            return False, "Thioester carbon connectivity is invalid"
    
        # Find the carbon adjacent to the carbonyl (alpha carbon)
        alpha_carbon_index = None
        for atom_index in neighboring_atoms:
          atom = mol.GetAtomWithIdx(atom_index)
          if atom.GetAtomicNum() == 6 and atom.GetTotalValence() == 4:
             alpha_carbon_index = atom_index
             break

        if alpha_carbon_index is None:
            return False, "Could not find the alpha carbon"
        
        alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_index)
        
        alpha_neighboring_atoms = [atom.GetIdx() for atom in alpha_carbon.GetNeighbors()]
        # find beta carbon
        beta_carbon_index = None
        for atom_index in alpha_neighboring_atoms:
            if atom_index != thioester_carbon_index:
                beta_carbon = mol.GetAtomWithIdx(atom_index)
                if beta_carbon.GetAtomicNum() == 6 and beta_carbon.GetTotalValence() == 3:
                  beta_carbon_index = atom_index
                  break

        if beta_carbon_index is None:
            return False, "Could not find beta carbon"
        
        # Find the bond between alpha and beta carbons
        bond = mol.GetBondBetweenAtoms(alpha_carbon_index,beta_carbon_index)
        
        if bond is None:
             return False, "No bond found between alpha and beta carbons"

        # Check for trans configuration.
        if bond.GetStereo() != Chem.BondStereo.STEREOE:
           return False, "Double bond not trans (E) configuration"

        if not bond.GetBondType() == Chem.BondType.DOUBLE:
          return False, "Bond between alpha and beta carbons is not double"

        # Check for a chain of at least 3 carbons after the double bond
        carbon_chain_pattern = Chem.MolFromSmarts("[CX4]-[CX4]-[CX4]")
        
        
        # Create a new molecule with the fatty acid part for counting carbons
        edit_mol = Chem.RWMol(mol)
        # Remove CoA
        for atom_idx in mol.GetSubstructMatch(coa_pattern):
              if mol.GetAtomWithIdx(atom_idx).GetSymbol() != 'S':
                  edit_mol.RemoveAtom(atom_idx)

        
        # Find a carbon chain attached to the double bond
        for match in edit_mol.GetSubstructMatches(carbon_chain_pattern):
            # check if the chain is in the acyl chain
            if mol.GetBondBetweenAtoms(match[0],match[1]) and mol.GetBondBetweenAtoms(match[1],match[2]):

              
               found_alpha = False
               for carbon_idx in match:
                 carbon = mol.GetAtomWithIdx(carbon_idx)
                 for neighbor_atom_idx in carbon.GetNeighbors():
                    if neighbor_atom_idx == alpha_carbon_index:
                         found_alpha = True
                         break
                 if found_alpha:
                    break
              
               if found_alpha:
                   
                  break
        else:
              return False, "No carbon chain found"

        # Check molecular weight - generally > 700 Da
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        if mol_wt < 700:
           return False, "Molecular weight too low for trans-2-enoyl-CoA"

    return True, "Molecule is a trans-2-enoyl-CoA"