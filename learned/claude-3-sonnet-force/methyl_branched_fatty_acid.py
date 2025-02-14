"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: CHEBI:35454 methyl-branched fatty acid
Any branched-chain fatty acid containing methyl branches only.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid is a branched-chain fatty acid containing methyl branches only.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for methyl branches
    methyl_pattern = Chem.MolFromSmarts("[CH3]")
    if not mol.HasSubstructMatch(methyl_pattern):
        return False, "No methyl groups found"

    # Check if all branches are methyl
    for atom in mol.GetAtoms():
        if atom.GetDegree() > 3 and atom.GetAtomicNum() == 6:
            subs = atom.GetHybridization()
            neighbors = [mol.GetAtomWithIdx(nei).GetAtomicNum() for nei in atom.GetNeighbors()]
            if subs == Chem.HybridizationType.SP3 and 6 in neighbors:
                if 1 not in neighbors:
                    return False, "Non-methyl branch found"

    # Check if molecule is acyclic
    if not Chem.Lipinski.IsAliphaticMolecule(mol):
        return False, "Molecule is not acyclic"

    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chain too short to be a fatty acid"

    return True, "Molecule is a methyl-branched fatty acid"