"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: CHEBI:18148 catecholamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine is a 4-(2-aminoethyl)pyrocatechol (4-(2-aminoethyl)benzene-1,2-diol)
    or a derivative formed by substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for pyrocatechol (benzene-1,2-diol) core
    pyrocatechol_pattern = Chem.MolFromSmarts("c1c(O)cc(O)cc1")
    pyrocatechol_matches = mol.GetSubstructMatches(pyrocatechol_pattern)
    if not pyrocatechol_matches:
        return False, "No pyrocatechol core found"

    # Look for 2-aminoethyl side chain
    ethylamine_pattern = Chem.MolFromSmarts("CCN")
    ethylamine_matches = mol.GetSubstructMatches(ethylamine_pattern)
    if not ethylamine_matches:
        return False, "No 2-aminoethyl side chain found"

    # Check connectivity between pyrocatechol and ethylamine
    connected = any(AllChem.GetShortestPath(mol, pyrocatechol_match[0], ethylamine_match[0])
                    for pyrocatechol_match in pyrocatechol_matches
                    for ethylamine_match in ethylamine_matches)
    if not connected:
        return False, "Pyrocatechol and ethylamine not connected"

    # Check for substitutions on the pyrocatechol core
    pyrocatechol_core_atom = mol.GetAtomWithIdx(pyrocatechol_matches[0][0])
    if pyrocatechol_core_atom.GetDegree() > 3:
        substituted = True
    else:
        substituted = False

    # Check for additional hydroxyl groups
    num_hydroxyl_groups = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
    if num_hydroxyl_groups > 2:
        substituted = True

    # Check molecular weight - catecholamines typically <300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 300:
        return False, "Molecular weight too high for catecholamine"

    if substituted:
        return True, "Substituted catecholamine derivative"
    else:
        return True, "Contains pyrocatechol core with 2-aminoethyl side chain"