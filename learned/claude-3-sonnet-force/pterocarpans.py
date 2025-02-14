"""
Classifies: CHEBI:26377 pterocarpans
"""
"""
Classifies: CHEBI:51811 pterocarpan
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    Pterocarpans are members of the class of benzofurochromene with a
    6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pterocarpan, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for pterocarpan scaffold pattern
    scaffold_pattern = Chem.MolFromSmarts("[C@]12Oc3ccccc3[C@@]1(O)COc1ccccc21")
    if mol.HasSubstructMatch(scaffold_pattern):
        # Check for fused rings
        fused_rings = AllChem.GetSymmSSSR(mol)
        if len(fused_rings) >= 4:
            return True, "Contains the pterocarpan scaffold with required fused rings"
        else:
            return False, "Does not have enough fused rings"
    else:
        # Check for substituted pterocarpan scaffold
        substituted_pattern = Chem.MolFromSmarts("[C@]12Oc3c(*)cccc3[C@@]1(O)COc1c(**)cccc21")
        if mol.HasSubstructMatch(substituted_pattern):
            fused_rings = AllChem.GetSymmSSSR(mol)
            if len(fused_rings) >= 4:
                return True, "Contains a substituted pterocarpan scaffold with required fused rings"
            else:
                return False, "Does not have enough fused rings"
        else:
            return False, "Does not contain the pterocarpan scaffold"