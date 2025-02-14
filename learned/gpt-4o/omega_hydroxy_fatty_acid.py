"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    An omega-hydroxy fatty acid has a carboxyl group at position 1 and a hydroxyl at the omega position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an omega-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    if not mol:
        return False, "Cannot parse the molecule from the SMILES"

    # Find carboxyl group
    carboxyl_group = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxyl_group):
        return False, "No carboxyl group (-COOH) found"

    # Identify the longest carbon chain containing the carboxyl group
    chain_mols = rdmolops.GetMolFrags(mol, asMols=True)
    longest_chain = max(chain_mols, key=lambda m: m.GetNumAtoms())

    # Ensure a single carboxyl group is at one end of the chain
    matches = longest_chain.GetSubstructMatches(carboxyl_group)
    if len(matches) != 1:
        return False, "More than one carboxyl group found"

    carboxyl_atom_index = matches[0][0]
    
    # Identify hydroxyl groups
    hydroxyl_group = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = longest_chain.GetSubstructMatches(hydroxyl_group)
    if not hydroxyl_matches:
        return False, "No hydroxyl groups found"

    # Check if the -OH is bonded to the terminal carbon farthest from the carboxyl group
    omega_hydroxyl_found = False
    
    for match in hydroxyl_matches:
        hydroxyl_index = match[0]
        hydroxyl_atom = longest_chain.GetAtomWithIdx(hydroxyl_index)
        for bond in hydroxyl_atom.GetBonds():
            bonded_atom = bond.GetOtherAtom(hydroxyl_atom)
            if bonded_atom.GetAtomicNum() == 6:  # Carbon
                terminal_carbon_index = bonded_atom.GetIdx()
                if terminal_carbon_index != carboxyl_atom_index:
                    omega_hydroxyl_found = True
                    break

    if not omega_hydroxyl_found:
        return False, "Hydroxyl group is not at the omega position of the longest chain"

    return True, "Contains a terminal carboxyl group and an omega-position hydroxyl group, fitting the omega-hydroxy fatty acid structure"