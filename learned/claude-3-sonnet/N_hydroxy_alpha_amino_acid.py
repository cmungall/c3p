"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: N-hydroxy-alpha-amino-acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino acid based on its SMILES string.
    These compounds have a carboxylic acid group and at least one hydrogen attached 
    to the amino group is replaced by a hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-hydroxy-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES and sanitize
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol is None:
            return False, "Invalid SMILES string"
    except:
        return False, "Invalid SMILES string"

    # Define key SMARTS patterns
    patterns = {
        # Alpha-amino acid with N-hydroxy group (more specific pattern)
        'n_hydroxy_aa': "[OX2H1,OX2H0]-[NX3](-[#1,*])-[CX4H1](-[#1,*])-C(=O)[OH1,O-]",
        # N,N-dihydroxy pattern
        'n_dihydroxy_aa': "[OX2H1,OX2H0]-[NX3](-[OX2H1,OX2H0])-[CX4H1](-[#1,*])-C(=O)[OH1,O-]",
        # Hydroxyimino patterns
        'hydroxyimino_aa': "[CX4H1](-[#1,*])(-[NX3](-[#1,*])-C(=N-[OX2H1,OX2H0])-[#1,*])-C(=O)[OH1,O-]",
        # Exclude patterns - functional groups that should not be present
        'ester': "[#6]-O-C(=O)-[#6]",
        'amide': "[#6]-N-C(=O)-[#6]",
        'anhydride': "[#6]-C(=O)-O-C(=O)-[#6]"
    }

    # Convert patterns to RDKit molecules
    queries = {name: Chem.MolFromSmarts(pattern) for name, pattern in patterns.items()}

    # Check for presence of core N-hydroxy amino acid structure
    is_n_hydroxy = mol.HasSubstructMatch(queries['n_hydroxy_aa'])
    is_n_dihydroxy = mol.HasSubstructMatch(queries['n_dihydroxy_aa'])
    is_hydroxyimino = mol.HasSubstructMatch(queries['hydroxyimino_aa'])

    if not (is_n_hydroxy or is_n_dihydroxy or is_hydroxyimino):
        return False, "No N-hydroxy-alpha-amino acid structure found"

    # Count the number of carboxylic acid groups
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH1,O-]")
    carboxyl_count = len(mol.GetSubstructMatches(carboxyl_pattern))
    
    if carboxyl_count == 0:
        return False, "No carboxylic acid group found"
    
    # Additional checks for false positives
    for exclude_pattern in ['ester', 'amide', 'anhydride']:
        if mol.HasSubstructMatch(queries[exclude_pattern]):
            matches = mol.GetSubstructMatches(queries[exclude_pattern])
            # Allow if the match doesn't involve the alpha-amino acid part
            if len(matches) > 1:
                return False, f"Contains {exclude_pattern} groups"

    # Determine specific type of N-hydroxy amino acid
    if is_n_dihydroxy:
        return True, "Contains N,N-dihydroxy group on alpha-amino acid"
    elif is_hydroxyimino:
        return True, "Contains hydroxyimino group on alpha-amino acid"
    else:
        return True, "Contains N-hydroxy group on alpha-amino acid"