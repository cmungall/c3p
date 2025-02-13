"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:36866 short-chain fatty acyl-CoA

A short-chain fatty acyl-CoA is defined as:
A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A
with the carboxy group of any short-chain fatty acid.

Short-chain fatty acids are typically defined as having 4-6 carbon atoms.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for coenzyme A substructure
    coenzyme_a_pattern = Chem.MolFromSmarts("CCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coenzyme_a_pattern):
        return False, "Coenzyme A substructure not found"

    # Look for acyl group attached to sulfur
    acyl_pattern = Chem.MolFromSmarts("[S](=O)(=O)CC")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl group attached to sulfur"

    # Count carbon atoms in acyl group
    acyl_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and atom.GetIsAromatic() and atom.GetTotalNumHs() == 0:
            continue  # Ignore aromatic carbons
        if atom.GetSymbol() == "C" and atom.GetIsInRingSize() > 6:
            continue  # Ignore large ring carbons
        if atom.GetSymbol() == "C" and atom.GetBonds()[0].GetBondType() == Chem.BondType.SINGLE:
            acyl_carbons.append(atom)

    n_acyl_carbons = len(acyl_carbons)
    if n_acyl_carbons < 4 or n_acyl_carbons > 6:
        return False, f"Acyl group has {n_acyl_carbons} carbons, should be between 4 and 6 for short-chain fatty acyl-CoA"

    return True, "Contains coenzyme A substructure with a short-chain fatty acyl group (4-6 carbons) attached to the sulfur"