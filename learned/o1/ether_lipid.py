"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: CHEBI:64675 ether lipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid is a lipid similar in structure to a glycerolipid but in which
    one or more of the carbon atoms on glycerol is bonded to an alkyl chain via an
    ether linkage, as opposed to the usual ester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the glycerol backbone pattern (3 connected carbons)
    glycerol_pattern = Chem.MolFromSmarts("[C;!$(C=[O,N,S])]-[C;!$(C=[O,N,S])]-[C;!$(C=[O,N,S])]")

    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone found"

    # Flag to check if at least one ether linkage is present
    has_ether_linkage = False

    for match in matches:
        glycerol_atoms = [mol.GetAtomWithIdx(idx) for idx in match]

        # Check if the three carbons are connected linearly
        if not (glycerol_atoms[0].IsInRing() or glycerol_atoms[1].IsInRing() or glycerol_atoms[2].IsInRing()):
            neighbors0 = set([a.GetIdx() for a in glycerol_atoms[0].GetNeighbors()])
            neighbors1 = set([a.GetIdx() for a in glycerol_atoms[1].GetNeighbors()])
            neighbors2 = set([a.GetIdx() for a in glycerol_atoms[2].GetNeighbors()])
            
            # Ensure linear connectivity
            if glycerol_atoms[1].GetIdx() in neighbors0 and glycerol_atoms[1].GetIdx() in neighbors2:
                # Check for ether linkages on glycerol carbons
                for atom in glycerol_atoms:
                    # Get oxygen neighbors
                    oxygen_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]

                    for oxygen in oxygen_neighbors:
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), oxygen.GetIdx())
                        if bond.GetBondType() == Chem.BondType.SINGLE:
                            # Check if oxygen is connected to an sp3 carbon (ether linkage)
                            other_carbon = None
                            for nbr in oxygen.GetNeighbors():
                                if nbr.GetIdx() != atom.GetIdx() and nbr.GetAtomicNum() == 6:
                                    # Ensure it's not part of a carbonyl group (ester linkage)
                                    # Check if this carbon is sp3 hybridized and not double bonded to oxygen
                                    nbr_bonds = nbr.GetBonds()
                                    is_carbonyl = any(bond.GetBondType() == Chem.BondType.DOUBLE and
                                                      bond.GetOtherAtom(nbr).GetAtomicNum() == 8
                                                      for bond in nbr_bonds)
                                    if not is_carbonyl:
                                        # Found ether linkage
                                        has_ether_linkage = True
                                        break
                        if has_ether_linkage:
                            break
                    if has_ether_linkage:
                        break
            if has_ether_linkage:
                break

    if not has_ether_linkage:
        return False, "No ether linkage found on glycerol backbone"

    return True, "Contains glycerol backbone with at least one alkyl chain attached via ether linkage"