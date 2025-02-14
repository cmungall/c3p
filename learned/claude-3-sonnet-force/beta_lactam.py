"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: CHEBI:35460 beta-lactam
A lactam in which the amide bond is contained within a four-membered ring,
which includes the amide nitrogen and the carbonyl carbon.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 4-membered ring with N and C=O groups
    beta_lactam_pattern = Chem.MolFromSmarts("[NR2][CR2]1[CR2][CR2]1=O")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Ensure the N and C=O are part of the same ring
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        n_atom = next((atom for atom in ring_atoms if atom.GetAtomicNum() == 7), None)
        c_atom = next((atom for atom in ring_atoms if atom.GetFormalCharge() == 0 and atom.GetAtomicNum() == 6 and atom.GetIsAromatic() == False), None)
        if n_atom and c_atom and any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in c_atom.GetBonds()):
            return True, "Contains a 4-membered ring with N and C=O groups"

    return False, "No valid beta-lactam ring found"