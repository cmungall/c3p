"""
Classifies: CHEBI:48953 cyclohexenones
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is defined as any six-membered alicyclic ketone having one double bond in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find six-membered rings
    ring_info = mol.GetRingInfo()
    for ring in ring_info.BondRings():
        if len(ring) == 6:
            # Check if it has exactly one ketone group (C=O)
            ketone_count = 0
            double_bond_count = 0

            for bond_idx in ring:
                bond = mol.GetBondWithIdx(bond_idx)
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    double_bond_count += 1
                    # Check if it's a double bond with oxygen (part of ketone)
                    begin_atom = bond.GetBeginAtom()
                    end_atom = bond.GetEndAtom()
                    if (begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 8) or \
                       (begin_atom.GetAtomicNum() == 8 and end_atom.GetAtomicNum() == 6):
                        ketone_count += 1

            if ketone_count == 1 and double_bond_count == 2:
                return True, "Contains a six-membered ring with one ketone group and another double bond"
    
    return False, "Does not meet cyclohexenone criteria"