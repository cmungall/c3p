"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
Here's a Python program to classify chemical entities as unsaturated fatty acyl-CoA:

"""
Classifies: CHEBI:35653 unsaturated fatty acyl-CoA

A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any unsaturated fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem, rdMolDescriptors

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for CoA backbone pattern
    coa_pattern = Chem.MolFromSmarts("OP(OP(OCC(C(=O)NCCC(=O)NCCSC(=O))O)=O)(OP(O)(OC[C@H]1O[C@@H]([C@@H](O)[C@H]1OP(O)(O)=O)N1C=NC2=C1N=CN=C2N)(O)=O)(O)=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA backbone found"

    # Look for unsaturated alkyl chain (at least one double bond)
    alkyl_chain = rdFMCS.FragmentMolWrap(mol, pattern=Chem.MolFromSmarts("[C;H3,H2,H1]-[C]-[C;H2]=[C;H2]-[C]"), completeRingsOnly=False)
    if not alkyl_chain:
        return False, "No unsaturated alkyl chain found"

    # Check for alkyl chain length (>= 6 carbons)
    alkyl_length = sum(1 for atom in alkyl_chain.mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if alkyl_length < 6:
        return False, "Alkyl chain too short for fatty acyl-CoA (< 6 carbons)"

    # Count double bonds in alkyl chain
    num_double_bonds = sum(1 for bond in alkyl_chain.mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if num_double_bonds == 0:
        return False, "No double bonds found in alkyl chain"

    # Check molecular weight (>= 500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for unsaturated fatty acyl-CoA (< 500 Da)"

    return True, "Contains CoA backbone and an unsaturated alkyl chain (>= 6 carbons)"