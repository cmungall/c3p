"""
Classifies: CHEBI:17761 ceramide
"""
"""
Classifies: CHEBI:36689 ceramide
A ceramide is an N-acyl-sphingoid base derivative with an amide-linked fatty acid. The fatty acids are typically saturated or monounsaturated with chain lengths from 14 to 26 carbon atoms; the presence of a hydroxyl group on carbon 2 is fairly common. Ceramides are generally precursors of more complex sphingolipids.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sphingosine backbone (C-N-C-C-C-OH)
    sphingosine_pattern = Chem.MolFromSmarts("[CH2X4][NX3][CH2X4][CH2X4][CH2X4][CH2X4][OH1]")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone found"
    
    # Look for amide-linked fatty acid (-C(=O)-N-)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide-linked fatty acid found"

    # Check fatty acid chain length (14-26 carbons)
    fatty_acid_atoms = mol.GetSubstructMatches(amide_pattern)[0]
    fatty_acid_chain = Chem.FindAllPathsOfLengthMult(mol, fatty_acid_atoms[0], 10, 26, useBondOrder=False)[0]
    if len(fatty_acid_chain) < 14 or len(fatty_acid_chain) > 26:
        return False, f"Fatty acid chain length ({len(fatty_acid_chain)}) not in range 14-26"

    # Check for hydroxyl group on C2 of fatty acid chain
    c2_atom = mol.GetAtomWithIdx(fatty_acid_chain[2])
    if c2_atom.GetTotalNumHs() != 1:
        return True, "Ceramide without hydroxyl group on C2 of fatty acid chain"

    return True, "Contains sphingosine backbone with amide-linked fatty acid chain of appropriate length and hydroxyl group on C2"