"""
Classifies: CHEBI:33913 corrinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    A corrinoid is a derivative of the corrin nucleus, which contains four reduced or partly
    reduced pyrrole rings joined in a macrocycle by three =C- groups and one direct
    carbon-carbon bond linking alpha positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
         bool: True if molecule is a corrinoid, False otherwise
         str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Pyrrole rings (reduced or partially reduced)
    # Note the [NX3;!$(N-)] which specifies a nitrogen with 3 neighbors, which cannot be in a N- bond, i.e. part of a ring. 
    # This indicates that it is part of a pyrrole not a pyridinium.
    # The corresponding carbons must have 3 neighbors - no double bonds to make it reduced or partially reduced
    pyrrole_pattern = Chem.MolFromSmarts("[NX3;!$(N-)]1[CX3H][CX3H]=[CX3H][CX3H]1") 
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    if len(pyrrole_matches) < 4:
         return False, f"Found {len(pyrrole_matches)} pyrrole rings. Need at least 4."
    
    # 2. Methine bridges (=CH-)
    # methine_pattern = Chem.MolFromSmarts("[CX3H]=[CX3H]([#6])[CX3H]=[CX3H]")
    methine_pattern = Chem.MolFromSmarts("[CX3H]=[CX3H][CX3H]=[CX3H]")
    methine_matches = mol.GetSubstructMatches(methine_pattern)
    if len(methine_matches) < 3:
        return False, f"Found {len(methine_matches)} methine bridges. Need at least 3."

    # 3. Direct C-C bond between alpha positions (adjacent to nitrogen in pyrrole)
    # this requires the definition of an alpha position on a pyrrole
    alpha_carbon_pattern = Chem.MolFromSmarts("[NX3;!$(N-)]1-[CX3H]1")
    alpha_matches = mol.GetSubstructMatches(alpha_carbon_pattern)
    if len(alpha_matches) < 4:
         return False, f"Found {len(alpha_matches)} alpha carbons on pyrrole rings. Need at least 4."
    
    # now we look for 2 alpha carbons connected
    direct_cc_pattern = Chem.MolFromSmarts("[CX3H]~[CX3H]")
    direct_cc_matches = mol.GetSubstructMatches(direct_cc_pattern)
    # check if any of the matches are connecting alpha carbons on pyrroles
    found_direct_cc = False
    for match in direct_cc_matches:
        for a1 in alpha_matches:
            for a2 in alpha_matches:
                if a1[0] in match and a2[0] in match and a1 != a2:
                    found_direct_cc = True
                    break
            if found_direct_cc:
                break
        if found_direct_cc:
            break


    if not found_direct_cc:
            return False, "No direct C-C bond between alpha carbons found."
   

    # If all criteria are met, it is likely a corrinoid
    return True, "Contains a corrin nucleus with four pyrroles, 3 methine bridges and a direct C-C bond."