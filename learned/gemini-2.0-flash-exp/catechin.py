"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    Catechins are flavan-3-ols, having the core structure of two benzene rings
    linked by a pyran ring and a hydroxyl group at the 3-position of the pyran ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a catechin, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core flavan-3-ol substructure: 2H-chromen with 3-OH and a phenyl group at position 2
    # '[C]1[C]([C]([O][C]2[C]1=[C]([C]=[C][C]=C2)[C]3=[C][C]=[C][C]=[C]3)[O])[H]' could be more specific, but also exclude other variants.
    flavan_3ol_pattern = Chem.MolFromSmarts("C1CC(O)C(Oc2cc(cc12)c3ccccc3)O")
    
    # Also try another pattern (more specific and allows for different bonds to OH group)
    #flavan_3ol_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](Cc2c(cc(cc2)O[C@H]1c1ccccc1)O)O")
    if not mol.HasSubstructMatch(flavan_3ol_pattern):
            return False, "Core flavan-3-ol structure not found."

    # Additional hydroxyl groups (allow multiple on each ring)
    # This can be too permissive
    #hydroxyphenyl_pattern = Chem.MolFromSmarts("c1cc(O)ccc1") # This pattern can match to any benzene ring, making this test too permissive
    #hydroxy_count = len(mol.GetSubstructMatches(hydroxyphenyl_pattern)) # count the phenyls
    #if hydroxy_count < 1:  # Allow at least one hydroxyl group on the phenyl ring. We already captured one in the flavan-3-ol_pattern
    #    return False, f"Too few hydroxyl groups. Found only {hydroxy_count}"

    # Acceptable substitutions. The main ring should have at least 2 OH groups, in addition to the 3-OH
    # Additional rings are also common, so it is hard to put a limit.
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, "Too few oxygen atoms"

    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, "Too few carbon atoms"
    

    return True, "Molecule contains the core flavan-3-ol structure with additional hydroxyl groups and possible substitutions."