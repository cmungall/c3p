"""
Classifies: CHEBI:24279 glucosinolate
"""
"""
Classifies: glucosinolate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdqueries

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    A glucosinolate is defined as 'Water-soluble anionic substituted thioglucosides.
    Glucosinolates have a central C atom which is bonded via an S atom to a glycone group
    and via an N atom to a sulfonated oxime group, and which also carries a side-group.
    The side-chain and sulfate group have an anti stereochemical configuration across
    the C=N double bond.'

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosinolate, False otherwise
        str: Reason for classification
    """

    from rdkit import Chem
    from rdkit.Chem import AllChem

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the thioglucose moiety connected via sulfur to a carbon
    # Thioglucose SMARTS pattern: a glucose ring connected via S to any group
    thioglucose_pattern = Chem.MolFromSmarts('O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O1)CS')
    if not mol.HasSubstructMatch(thioglucose_pattern):
        return False, "No thioglucose moiety connected via sulfur found"

    # Look for the central carbon atom bonded to sulfur and nitrogen
    central_carbon_pattern = Chem.MolFromSmarts('[C;X4](S[*])[N]=[O]')
    central_carbon_matches = mol.GetSubstructMatches(central_carbon_pattern)
    if len(central_carbon_matches) == 0:
        return False, "No central carbon bonded to sulfur and nitrogen found"

    # Look for sulfonated oxime group attached to nitrogen
    sulfonated_oxime_pattern = Chem.MolFromSmarts('N=O-S(=O)(=O)[O-]')
    if not mol.HasSubstructMatch(sulfonated_oxime_pattern):
        return False, "No sulfonated oxime group found"

    # Check that the sulfur of thioglucose is connected to the central carbon
    # and the nitrogen of the sulfonated oxime is connected to the same carbon
    # We can check the bonding pattern around the central carbon
    is_glucosinolate = False
    for match in central_carbon_matches:
        c_idx = match[0]  # Index of the central carbon
        atom = mol.GetAtomWithIdx(c_idx)
        neighbors = [nbr.GetAtomicNum() for nbr in atom.GetNeighbors()]
        if neighbors.count(16) >=1 and neighbors.count(7) >=1:
            is_glucosinolate = True
            break
    if not is_glucosinolate:
        return False, "Central carbon is not bonded to both sulfur and nitrogen"

    # Check for the presence of a side-chain (R group)
    # The central carbon should have four bonds: S, N, R-group, and possibly a hydrogen
    atom = mol.GetAtomWithIdx(c_idx)
    num_neighbors = len(atom.GetNeighbors())
    if num_neighbors < 3:
        return False, "Central carbon does not have enough substituents for a side-chain"

    # Optionally, verify the C=N double bond
    # Check if the bond between the central carbon and nitrogen is a double bond
    bond = mol.GetBondBetweenAtoms(c_idx, match[2])  # match[2] is the nitrogen
    if bond is None or bond.GetBondType() != Chem.BondType.DOUBLE:
        return False, "No C=N double bond found between central carbon and nitrogen"

    return True, "Molecule matches the structural features of a glucosinolate"


__metadata__ = {
    'chemical_class': {
        'name': 'glucosinolate',
        'definition': 'Water-soluble anionic substituted thioglucosides. Glucosinolates have a central C atom which is bonded via an S atom to a glycone group and via an N atom to a sulfonated oxime group, and which also carries a side-group. The side-chain and sulfate group have an anti stereochemical configuration across the C=N double bond.'
    },
    'references': ['CHEBI:24278'],
}