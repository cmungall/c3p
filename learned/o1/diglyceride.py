"""
Classifies: CHEBI:18035 diglyceride
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride based on its SMILES string.
    A diglyceride is a glycerol backbone with two fatty acid chains attached via ester bonds,
    and the third position can be either a hydroxyl group (-OH) or an alkyl ether.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diglyceride, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens to the molecule for accurate matching
    mol = Chem.AddHs(mol)

    # Define glycerol backbone pattern with three carbons each attached to oxygen
    glycerol_pattern = Chem.MolFromSmarts("[C;H2][C;H][C;H2]")
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone found"

    # Iterate over matches to find the correct glycerol backbone
    for match in matches:
        c1_idx, c2_idx, c3_idx = match[0], match[1], match[2]

        # Get oxygen atoms attached to each carbon
        o1 = None
        o2 = None
        o3 = None

        for nbr in mol.GetAtomWithIdx(c1_idx).GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                o1 = nbr
        for nbr in mol.GetAtomWithIdx(c2_idx).GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                o2 = nbr
        for nbr in mol.GetAtomWithIdx(c3_idx).GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                o3 = nbr

        if not o1 or not o2 or not o3:
            continue  # Incomplete glycerol backbone

        # Check esterification status of oxygen atoms
        ester_count = 0
        free_oxygen_count = 0

        for o_atom in [o1, o2, o3]:
            o_idx = o_atom.GetIdx()
            bonds = o_atom.GetBonds()
            is_ester = False
            is_free = False

            for bond in bonds:
                bond_type = bond.GetBondType()
                other_atom = bond.GetOtherAtom(o_atom)

                if bond_type == Chem.BondType.SINGLE and other_atom.GetAtomicNum() == 6:
                    # Check if connected to a carbonyl carbon (ester linkage)
                    for obond in other_atom.GetBonds():
                        if obond.GetBondType() == Chem.BondType.DOUBLE and obond.GetOtherAtom(other_atom).GetAtomicNum() == 8:
                            is_ester = True
                            break
                elif bond_type == Chem.BondType.SINGLE and other_atom.GetAtomicNum() == 1:
                    # Connected to hydrogen (hydroxyl group)
                    is_free = True
                elif bond_type == Chem.BondType.SINGLE and other_atom.GetAtomicNum() == 6:
                    # Could be ether linkage
                    is_free = True

                if is_ester:
                    break

            if is_ester:
                ester_count += 1
            elif is_free:
                free_oxygen_count += 1
            else:
                # Oxygen is neither esterified nor free hydroxyl/ether
                break

        if ester_count == 2 and free_oxygen_count == 1:
            # Optionally, check that fatty acid chains are long enough
            # This can be done by checking the chain lengths attached to esterified oxygens
            return True, "Contains glycerol backbone with two fatty acid chains attached via ester bonds and one hydroxyl or ether group"

    return False, "Does not match diglyceride structure"