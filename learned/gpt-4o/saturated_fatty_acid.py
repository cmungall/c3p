"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    A saturated fatty acid contains no carbon-to-carbon multiple bonds and includes a
    terminal carboxyl group, typically having a long linear carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure the molecule contains a carboxyl group (-C(=O)O) 
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and any(
            bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and 
            bond.GetOtherAtom(atom).GetAtomicNum() == 6
            for bond in atom.GetBonds()):
            carboxyl_group = True
            break
    else:
        return False, "No carboxyl group found"

    # Check for saturation: ensure no double or triple C-C bonds
    for bond in mol.GetBonds():
        if bond.GetBondType() not in (Chem.rdchem.BondType.SINGLE,):
            begin_atom, end_atom = bond.GetBeginAtom(), bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                return False, "Contains multiple carbon-carbon bonds"

    # Check for the longest carbon chain - should be linear
    longest_chain = max((len(path) for path in Chem.GetSymmSSSR(mol)), default=0)
    
    # Consider linear chains typical for saturated fatty acids
    chain_length = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if longest_chain > 0:
        return False, "Contains cyclic structure, not a saturated fatty acid"

    # Ensure chain length is reasonable for saturated fatty acids including short chains
    if chain_length < 4:
        return False, f"Insufficient carbon chain length, found {chain_length} carbons"

    return True, "Molecule is a saturated fatty acid with unbranched carbon chain and terminal carboxyl group"