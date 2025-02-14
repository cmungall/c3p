"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
"""
Classifies: 11,12-saturated fatty acyl-CoA(4-)
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is an 11,12-saturated fatty acyl-CoA(4-) based on its SMILES string.
    An 11,12-saturated fatty acyl-CoA(4-) is a fatty acyl-CoA(4-) molecule in which the bond between
    carbons 11 and 12 in the fatty acyl chain is saturated (a single bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11,12-saturated fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure molecule has explicit hydrogens
    mol = Chem.AddHs(mol)

    # Check for CoA substructure (simplified pattern)
    coa_smarts = 'NC(=O)CCNC(=O)[C@H](O)C(C)(C)CO[P](=O)([O-])O[P](=O)([O-])OC[C@H]1O[C@H](C)[C@H](O)[C@@H]1O'
    coa_mol = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_mol):
        return False, "Coenzyme A moiety not found"

    # Find the thioester linkage (C(=O)S)
    thioester_smarts = 'C(=O)S'
    thioester_mol = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatch(thioester_mol)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Identify the fatty acyl chain starting from the carbonyl carbon
    fatty_acyl_start = thioester_matches[0]  # carbonyl carbon index

    # Traverse the fatty acyl chain to number the carbons
    chain_atoms = [fatty_acyl_start]
    current_atom = mol.GetAtomWithIdx(fatty_acyl_start)
    prev_atom = None
    while len(chain_atoms) < 20:  # Limit to avoid infinite loops
        neighbors = [nbr for nbr in current_atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev_atom]
        if not neighbors:
            break
        next_atom = neighbors[0]
        chain_atoms.append(next_atom.GetIdx())
        prev_atom = current_atom.GetIdx()
        current_atom = next_atom

    if len(chain_atoms) < 12:
        return False, "Fatty acyl chain is shorter than 12 carbons"

    # Check the bond between carbons 11 and 12
    c11_idx = chain_atoms[10]  # zero-based indexing
    c12_idx = chain_atoms[11]
    bond = mol.GetBondBetweenAtoms(c11_idx, c12_idx)
    if bond is None:
        return False, "No bond between carbons 11 and 12"
    bond_type = bond.GetBondType()
    if bond_type == Chem.rdchem.BondType.SINGLE:
        return True, "Bond between carbons 11 and 12 is saturated (single bond)"
    elif bond_type == Chem.rdchem.BondType.DOUBLE:
        return False, "Bond between carbons 11 and 12 is unsaturated (double bond)"
    else:
        return False, f"Bond between carbons 11 and 12 is neither single nor double bond (bond type: {bond_type})"