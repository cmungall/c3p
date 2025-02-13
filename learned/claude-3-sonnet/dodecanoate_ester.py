"""
Classifies: CHEBI:87659 dodecanoate ester
"""
"""
Classifies: CHEBI:32847 dodecanoate ester
Defined as: Any fatty acid ester in which the carboxylic acid component is lauric acid (dodecanoic acid).
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find dodecanoate fragment
    dodecanoate_fragment = Chem.MolFromSmiles("CCCCCCCCCCCC(=O)O")
    dodecanoate_matches = mol.GetSubstructMatches(dodecanoate_fragment)
    if not dodecanoate_matches:
        return False, "Dodecanoate fragment not found"

    # Check for ester bonds involving the dodecanoate fragment
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    dodecanoate_ester_matches = []
    for match in ester_matches:
        ester_atom = mol.GetAtomWithIdx(match[1])
        for dodecanoate_match in dodecanoate_matches:
            if ester_atom.GetIdx() in dodecanoate_match:
                dodecanoate_ester_matches.append(match)
                break

    if not dodecanoate_ester_matches:
        return False, "No ester bond involving the dodecanoate fragment found"

    # Additional checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for dodecanoate ester"

    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 15 or o_count < 2:
        return False, "Insufficient carbon and oxygen atoms for dodecanoate ester"

    return True, "Contains dodecanoate fragment as the carboxylic acid component of an ester bond"