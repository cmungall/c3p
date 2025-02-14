"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: trans-2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    A trans-2-enoyl-CoA is an unsaturated fatty acyl-CoA that results from the formal condensation
    of the thiol group of coenzyme A with the carboxy group of any 2,3-trans-enoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trans-2-enoyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define CoA pattern (simplified)
    coa_pattern_smiles = "NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@@H]1O[C@H](n2cnc3c(N)ncnc32)[C@H](O)[C@H]1OP(=O)(O)O"
    coa_pattern = Chem.MolFromSmiles(coa_pattern_smiles)
    if coa_pattern is None:
        return False, "Failed to generate CoA substructure"

    # Check for CoA substructure
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Define thioester linkage pattern: C(=O)-S
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Iterate over thioester matches
    for thioester_match in thioester_matches:
        # Indices in thioester_match correspond to [carbonyl carbon, sulfur atom]
        carbonyl_c_idx = thioester_match[0]
        sulfur_idx = thioester_match[1]

        # Check that sulfur is part of CoA by ensuring sulfur atom is included in the CoA match
        # Get the CoA substructure match
        coa_match = mol.GetSubstructMatch(coa_pattern)
        if sulfur_idx not in coa_match:
            continue  # Sulfur is not part of CoA, skip this thioester

        # Extract the acyl chain attached to the carbonyl carbon
        # Break bond between carbonyl carbon and sulfur atom to isolate acyl chain
        emol = Chem.RWMol(mol)
        emol.RemoveBond(carbonyl_c_idx, sulfur_idx)
        mol_frag = emol.GetMol()

        # Get the fragments after bond removal
        frags = Chem.GetMolFrags(mol_frag, asMols=True, sanitizeFrags=False)

        # Identify the acyl chain fragment (containing the carbonyl carbon)
        acyl_chain = None
        for frag in frags:
            if frag.HasSubstructMatch(Chem.MolFromSmiles("C=O")):
                acyl_chain = frag
                break

        if acyl_chain is None:
            continue  # Could not find acyl chain fragment

        # Check for trans double bond between C2 and C3 in the acyl chain
        # Define a pattern for the acyl chain starting with carbonyl carbon
        # and having a trans double bond between C2 and C3
        trans_double_bond_pattern = Chem.MolFromSmarts("C(=O)C\C=C\C")

        if acyl_chain.HasSubstructMatch(trans_double_bond_pattern):
            return True, "Contains CoA moiety with acyl chain having trans double bond between C2 and C3"

    return False, "Does not contain acyl chain with trans double bond between C2 and C3 attached to CoA"

__metadata__ = {
    'chemical_class': {
        'id': '',
        'name': 'trans-2-enoyl-CoA',
        'definition': 'An unsaturated fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any 2,3-trans-enoic acid.',
        'parents': []
    },
    'message': None,
    'success': True,
    'error': '',
}