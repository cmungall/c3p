"""
Classifies: CHEBI:17334 penicillin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    A penicillin is defined as having a penam core, two methyl substituents at position 2,
    a carboxylate substituent at position 3, and a carboxamido group at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Penam core pattern (with numbering for easy referencing)
    #     7    6     5
    #   --N--[C@H]--[CH]--
    #  |   |     |     |
    #  |   C------C------S
    #  |   2     3     4
    #   --    C=O  --
    #
    penam_core_pattern = Chem.MolFromSmarts("[N1][C@H]([C@H]2[S][C@](C)(C)[C@H]2C(O)=O)C(=O)1")

    if not mol.HasSubstructMatch(penam_core_pattern):
        return False, "Penam core not found"


    # 2. Two methyl substituents at position 2. These are already checked in the penam core definition above.
    # 3. Carboxylate substituent at position 3. This is already checked in the penam core definition above.

    # 4. Carboxamido group at position 6 (R-CONH-). It's not a specific R group.
    carboxamido_pattern = Chem.MolFromSmarts("N[C@H](C=O)")
    carboxamido_matches = mol.GetSubstructMatches(carboxamido_pattern)

    # Use the core substructure to count the matches, instead of the whole molecule.
    substructure_matches = mol.GetSubstructMatches(penam_core_pattern)
    if len(substructure_matches) == 0:
        return False, "Penam core not found when checking amido group."
    
    match = substructure_matches[0] # Get the atoms from the first match.
    
    sub_mol = Chem.PathToSubmol(mol,match)
    sub_match = sub_mol.GetSubstructMatches(carboxamido_pattern)
    
    if len(sub_match) !=1:
        return False, f"Carboxamido group at position 6 not found within the penam core. Found {len(sub_match)} matches."


    return True, "Molecule contains penam core, two methyl substituents, carboxylate, and carboxamido group"