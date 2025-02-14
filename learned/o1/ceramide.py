"""
Classifies: CHEBI:17761 ceramide
"""
"""
Classifies: ceramide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    A ceramide is an N-acyl-sphingoid base with an amide-linked fatty acid.
    The fatty acid is typically saturated or monounsaturated with 14 to 26 carbons.
    The sphingoid base is a long-chain amino alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ceramide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for amide bond (C(=O)-N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond found"

    # Assume first amide bond is the linkage
    amide_match = amide_matches[0]
    carbonyl_c_idx = amide_match[0]
    nitrogen_idx = amide_match[2]

    # Break the molecule at the amide bond to separate fatty acid and sphingoid base
    amide_bond = mol.GetBondBetweenAtoms(carbonyl_c_idx, nitrogen_idx)
    bond_idx = amide_bond.GetIdx()
    fragmented_mol = Chem.FragmentOnBonds(mol, [bond_idx])

    # Get the fragments
    frags = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)
    if len(frags) != 2:
        return False, "Could not separate fatty acid and sphingoid base"

    # Identify fatty acid and sphingoid base fragments
    fatty_acid = None
    sphingoid_base = None
    for frag in frags:
        atom_nums = [atom.GetAtomicNum() for atom in frag.GetAtoms()]
        num_carbons = atom_nums.count(6)
        num_nitrogens = atom_nums.count(7)
        num_oxygens = atom_nums.count(8)
        num_hydroxyls = len(frag.GetSubstructMatches(Chem.MolFromSmarts("[OX2H]")))

        if num_nitrogens == 0 and num_carbons >= 14:
            fatty_acid = frag
            fatty_acid_carbons = num_carbons
        elif num_nitrogens >= 1:
            sphingoid_base = frag
            sphingoid_base_carbons = num_carbons
            sphingoid_base_hydroxyls = num_hydroxyls

    if fatty_acid is None or sphingoid_base is None:
        return False, "Could not identify fatty acid and sphingoid base fragments"

    # Check fatty acid chain length
    if fatty_acid_carbons < 14 or fatty_acid_carbons > 26:
        return False, f"Fatty acid chain length is {fatty_acid_carbons}, which is outside 14-26 carbons"

    # Check if fatty acid is saturated or monounsaturated
    fatty_acid_unsaturations = rdMolDescriptors.CalcNumDoubleBonds(fatty_acid)
    if fatty_acid_unsaturations > 1:
        return False, f"Fatty acid has {fatty_acid_unsaturations} double bonds, should be saturated or monounsaturated"

    # Check sphingoid base chain length
    if sphingoid_base_carbons < 12:
        return False, f"Sphingoid base chain length is {sphingoid_base_carbons}, which is too short"

    # Check for hydroxyl groups in sphingoid base
    if sphingoid_base_hydroxyls < 1:
        return False, "No hydroxyl groups found on sphingoid base"
    else:
        hydroxyl_positions = sphingoid_base.GetSubstructMatches(Chem.MolFromSmarts("[C][O][H]"))
        if not hydroxyl_positions:
            return False, "Hydroxyl groups not attached to carbon in sphingoid base"

    # Optional: Check for hydroxyl group on carbon 2 (common but not mandatory)
    # This would require mapping atom indices which can be complex
    # Skipping for simplicity

    return True, "Molecule is a ceramide with appropriate chain lengths and functional groups"