"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids are endogenous ligands that activate cannabinoid receptors,
    typically consisting of long-chain polyunsaturated fatty acids connected
    via an amide linkage to ethanolamine or an ester/ether linkage to glycerol.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is an endocannabinoid, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Count the number of carbons in the molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 16:
        return False, "Molecule does not have a long carbon chain (less than 16 carbons)."

    # Count the number of double bonds
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bond_count +=1
    if double_bond_count < 2:
        return False, "Insufficient number of double bonds for polyunsaturation (less than 2)."

    # Define SMARTS patterns for functional groups
    # Amide linkage to ethanolamine: -C(=O)-N-CCO
    ethanolamine_amide_pattern = Chem.MolFromSmarts("C(=O)NCCO")
    # Ester linkage to glycerol: -C(=O)OCC(O)CO
    glycerol_ester_pattern = Chem.MolFromSmarts("C(=O)OCC(O)CO")
    # Ether linkage to glycerol: -COCC(O)CO
    glycerol_ether_pattern = Chem.MolFromSmarts("COCC(O)CO")

    # Check for amide linkage to ethanolamine
    if mol.HasSubstructMatch(ethanolamine_amide_pattern):
        return True, "Contains amide linkage to ethanolamine and polyunsaturated fatty acid chain."

    # Check for ester linkage to glycerol
    if mol.HasSubstructMatch(glycerol_ester_pattern):
        return True, "Contains ester linkage to glycerol and polyunsaturated fatty acid chain."

    # Check for ether linkage to glycerol
    if mol.HasSubstructMatch(glycerol_ether_pattern):
        return True, "Contains ether linkage to glycerol and polyunsaturated fatty acid chain."

    return False, "Molecule does not match common endocannabinoid structures."