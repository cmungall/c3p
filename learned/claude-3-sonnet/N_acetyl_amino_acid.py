"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: N-acetyl-amino acid
Definition: An N-acyl-amino acid that has acetyl as the acyl group
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an N-acetyl-amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (including deprotonated form)
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[OH,O-]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Look for N-acetyl group (CH3-C(=O)-N), allowing for hydrogen on N
    n_acetyl = Chem.MolFromSmarts('[CH3]C(=O)[NH1,NH0]')
    if not mol.HasSubstructMatch(n_acetyl):
        return False, "No N-acetyl group found"

    # Pattern for amino acid core with N-acetyl
    # Matches both linear and cyclic cases
    core_patterns = [
        # Linear amino acids: acetyl-N-C-C(=O)OH with any substitution on alpha carbon
        Chem.MolFromSmarts('[CH3]C(=O)[NH1,NH0]C([*,H])[*]C(=O)[OH,O-]'),
        # Cyclic amino acids (like proline derivatives)
        Chem.MolFromSmarts('[CH3]C(=O)N1[CH2,CH]CC[CH2,CH]1C(=O)[OH,O-]'),
        # Alternative connection pattern
        Chem.MolFromSmarts('[CH3]C(=O)[NH1,NH0]C[*]C(=O)[OH,O-]')
    ]
    
    core_found = False
    for pattern in core_patterns:
        if mol.HasSubstructMatch(pattern):
            core_found = True
            break
            
    if not core_found:
        return False, "No N-acetyl amino acid core structure found"

    # Verify connectivity between N-acetyl and carboxylic acid
    n_acetyl_matches = mol.GetSubstructMatches(n_acetyl)
    carboxyl_matches = mol.GetSubstructMatches(carboxylic_acid)
    
    for n_match in n_acetyl_matches:
        n_idx = n_match[2]  # Index of nitrogen atom in N-acetyl match
        for c_match in carboxyl_matches:
            c_idx = c_match[0]  # Index of carbon atom in carboxyl match
            path = Chem.GetShortestPath(mol, n_idx, c_idx)
            if path and 2 <= len(path) <= 8:  # Increased maximum path length
                # Additional check: Count atoms between N and COOH
                atoms_between = len([x for x in path[1:-1] 
                                  if mol.GetAtomWithIdx(x).GetAtomicNum() != 1])
                if 1 <= atoms_between <= 6:  # At least one carbon between N and COOH
                    return True, "Contains N-acetyl group properly connected to amino acid core"
                
    return False, "N-acetyl group not properly connected to amino acid core"