"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: isothiocyanate compounds (R-N=C=S)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    Isothiocyanates have the general formula R-N=C=S where R is an organic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isothiocyanate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for N=C=S group
    isothiocyanate_pattern = Chem.MolFromSmarts("[NX2]=[CX2]=[SX1]")
    matches = mol.GetSubstructMatches(isothiocyanate_pattern)
    
    if not matches:
        return False, "No isothiocyanate group found"
    
    # Get basic molecular properties
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    num_atoms = mol.GetNumAtoms()
    
    valid_groups = 0
    for match in matches:
        n_atom, c_atom, s_atom = [mol.GetAtomWithIdx(idx) for idx in match]
        
        # Validate nitrogen
        if n_atom.GetDegree() != 2 or n_atom.GetTotalValence() != 3:
            continue
            
        # Validate carbon
        if c_atom.GetDegree() != 2 or c_atom.GetTotalValence() != 4:
            continue
            
        # Validate sulfur
        if s_atom.GetDegree() != 1 or s_atom.GetTotalValence() != 2:
            continue
            
        # Get R group (atom connected to nitrogen that's not the central carbon)
        r_atoms = [atom for atom in n_atom.GetNeighbors() if atom.GetIdx() != c_atom.GetIdx()]
        if not r_atoms:
            continue
        r_atom = r_atoms[0]
        
        # Validate R group
        if r_atom.GetAtomicNum() not in [6, 7, 8]:  # Must be C, N, or O
            continue
            
        # Check bond types
        n_c_bond = mol.GetBondBetweenAtoms(n_atom.GetIdx(), c_atom.GetIdx())
        c_s_bond = mol.GetBondBetweenAtoms(c_atom.GetIdx(), s_atom.GetIdx())
        
        if not (n_c_bond.GetBondType() == Chem.BondType.DOUBLE and 
                c_s_bond.GetBondType() == Chem.BondType.DOUBLE):
            continue
        
        # Check if the isothiocyanate group is a significant feature
        # by comparing its size to the overall molecule
        ncs_size = 3  # N=C=S group size
        if ncs_size / num_atoms < 0.1:  # If group is less than 10% of molecule
            if mol_weight > 500:  # And molecule is large
                continue  # Likely a complex molecule where -NCS is not the main feature
        
        valid_groups += 1
    
    if valid_groups == 0:
        return False, "No valid isothiocyanate groups found"
        
    # Additional checks for problematic features
    problematic_patterns = [
        Chem.MolFromSmarts("[N+]"), # Quaternary nitrogen
        Chem.MolFromSmarts("[N-]"), # Negatively charged nitrogen
        Chem.MolFromSmarts("[S+]"), # Positively charged sulfur
        Chem.MolFromSmarts("[S-]"), # Negatively charged sulfur
    ]
    
    for pattern in problematic_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return False, "Contains charged species that indicate different chemical class"
    
    # Check if molecule is too complex
    if mol_weight > 800:
        return False, "Molecule too complex to be classified primarily as isothiocyanate"
    
    if valid_groups == 1:
        return True, "Contains one isothiocyanate group (R-N=C=S)"
    else:
        return True, f"Contains {valid_groups} isothiocyanate groups"