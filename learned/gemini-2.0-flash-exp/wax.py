"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    Waxes are generally long-chain esters or alcohols with flexibility
    and moderate to high molecular weight.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 0. Check for heteroatoms other than C, H, O
    heteroatoms = [atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6, 8]]
    if len(heteroatoms) > 0:
        return False, "Contains heteroatoms other than C, H and O; not typical for wax"

    # 1. Check for rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains rings; not typical for a wax"

    # 2. Check for long continuous carbon chains (at least 10)
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    long_chain_carbons = set()
    
    all_chains = []
    for start_atom in carbon_atoms:
      for end_atom in carbon_atoms:
          if start_atom == end_atom:
             continue
          path = Chem.GetShortestPath(mol, start_atom.GetIdx(), end_atom.GetIdx())
          if path is not None:
              is_all_carbon = all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in path)
              if is_all_carbon:
                  all_chains.append(path)
                  if len(path) >= 10:
                        for idx in path:
                            long_chain_carbons.add(idx)

    # Count total number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # 3. Check that a good percentage of the carbons is in long chains
    if len(long_chain_carbons) < 0.6 * c_count: # at least 60% should be in long chains
        return False, f"Too few carbons in long chains ({len(long_chain_carbons)} out of {c_count})."

    # 4. Check for ester linkages (or long-chain alcohols)
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4]") # Ester next to a carbon
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H1]") # long-chain alcohols check
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)

    if len(ester_matches) == 0 and len(alcohol_matches) == 0:
        return False, f"No long-chain ester or alcohol groups detected."

    # Check that the ester/alcohol group is part of the long chain. If it's an ester it must have another alkyl
    if len(ester_matches) > 0 :
        ester_match_ok = False
        for match in ester_matches:           
           ester_carbon1 = mol.GetAtomWithIdx(match[0])
           ester_carbon2 = mol.GetAtomWithIdx(match[3]) # carbon to which the ester oxygen is linked

           found_chain1=False
           found_chain2=False
           for chain in all_chains:
                if ester_carbon1.GetIdx() in chain:
                   found_chain1 = True
                if ester_carbon2.GetIdx() in chain:
                   found_chain2 = True

           if found_chain1 and found_chain2:
               ester_match_ok = True
               break
        if not ester_match_ok:
             return False, "Ester group is not part of long carbon chain"
            
    if len(alcohol_matches) > 0:
        alcohol_match_ok = False
        for match in alcohol_matches:
            alcohol_carbon = mol.GetAtomWithIdx(match[0])
            
            found_chain = False
            for chain in all_chains:
                if alcohol_carbon.GetIdx() in chain:
                    found_chain = True
                    break
            if found_chain:
                alcohol_match_ok = True
                break
        if not alcohol_match_ok:
             return False, "Alcohol group is not part of long carbon chain"


    # 5 Correct check number of oxygens (max 3)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    required_o_for_esters = 0
    if len(ester_matches) > 0:
       required_o_for_esters = len(ester_matches)*2
    required_o_for_alcohols = 0
    if len(alcohol_matches) > 0:
      required_o_for_alcohols = len(alcohol_matches)

    if o_count != required_o_for_esters + required_o_for_alcohols:
         return False, f"Inconsistent number of oxygens {o_count}, expected {required_o_for_esters+required_o_for_alcohols} based on ester/alcohol count."

    if o_count > 3:
        return False, f"Too many oxygen atoms ({o_count}), waxes should have at most 3."

    # 6. Check for rotatable bonds (at least 8, proportional to the number of carbon atoms)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 0.3 * c_count or n_rotatable < 8:
          return False, f"Too few rotatable bonds ({n_rotatable}). Should be at least 8 and about 30% of the carbon count."

    # 7. Check molecular weight (approx range 250-1200)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 1200:
        return False, f"Molecular weight out of range ({mol_wt}). Should be between 250 and 1200"

    # 8. Check if chains are mostly saturated
    unsaturated_bond_pattern = Chem.MolFromSmarts("C=C")
    unsaturated_bonds = mol.GetSubstructMatches(unsaturated_bond_pattern)
    
    if len(unsaturated_bonds) > c_count/10: # max 1 unsaturation per 10 carbons
         return False, "Too many unsaturations"

    return True, "Meets the criteria for a wax (long, largely saturated chains, ester or alcohol, flexible, good MW)."