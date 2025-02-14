"""
Classifies: CHEBI:75874 para-terphenyl
"""
"""
Classifies: CHEBI:39158 para-terphenyl
A ring assembly based on a 1,4-diphenylbenzene skeleton and its substituted derivatives thereof.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl based on its SMILES string.
    A para-terphenyl is a ring assembly based on a 1,4-diphenylbenzene skeleton
    and its substituted derivatives thereof.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a para-terphenyl, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for para-terphenyl core structure
    # C1=CC(=CC=C1C2=CC=CC=C2)C3=CC=CC=C3
    terphenyl_pattern = Chem.MolFromSmarts("c1ccc(-c2ccc(-c3ccccc3)cc2)cc1")
    if not mol.HasSubstructMatch(terphenyl_pattern):
        return False, "No para-terphenyl core structure found"

    # Check for allowed atoms
    allowed_atoms = [6, 8, 16, 17, 35, 53]  # C, O, S, Cl, Br, I
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Atom {atom.GetSymbol()} not allowed in para-terphenyls"

    # Check for common functional groups
    functional_groups = ["CO", "OC", "CO.O", "C(=O)O", "C(=O)C", "C=C"]
    for fg in functional_groups:
        fg_pattern = Chem.MolFromSmarts(fg)
        if mol.HasSubstructMatch(fg_pattern):
            break
    else:
        return False, "No allowed functional groups found"

    # Additional checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 1000:
        return False, "Molecular weight out of expected range for para-terphenyls"

    n_rings = mol.GetRingInfo().NumRings()
    if n_rings < 3 or n_rings > 6:
        return False, "Number of rings out of expected range for para-terphenyls"

    return True, "Meets structural and functional group requirements for para-terphenyls"