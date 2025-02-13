"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
"""
Classifies: N-acylsphinganine
A ceramide consisting of sphinganine in which one of the amino hydrogens 
is substituted by a fatty acyl group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an N-acylsphinganine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core structural requirements:
    # 1. Strict sphinganine core with correct stereochemistry
    # [C@H] and [C@@H] enforce specific stereochemistry
    # No double bonds in the core (distinguishes from sphingosine)
    sphinganine_core = Chem.MolFromSmarts("[CH2X4][C@H]([OX2H1])[C@H]([NX3H1])[CH2X4][OX2H1]")
    if not mol.HasSubstructMatch(sphinganine_core):
        return False, "Missing characteristic sphinganine core structure with correct stereochemistry"

    # 2. N-acyl group (amide) directly attached to the sphinganine core
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[#6]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    # Check if any amide is properly connected to the sphinganine nitrogen
    valid_amide = False
    for match in amide_matches:
        amide_c = mol.GetAtomWithIdx(match[0])
        for neighbor in amide_c.GetNeighbors():
            if neighbor.GetAtomicNum() == 7:  # Nitrogen
                valid_amide = True
                break
    if not valid_amide:
        return False, "No N-acyl group properly connected to sphinganine core"

    # 3. Verify no extra hydroxyl groups on the core
    # (distinguishes from phytosphingosine which has an extra OH)
    core_carbons = mol.GetSubstructMatches(sphinganine_core)[0]
    for c_idx in core_carbons:
        carbon = mol.GetAtomWithIdx(c_idx)
        oh_count = sum(1 for n in carbon.GetNeighbors() 
                      if n.GetAtomicNum() == 8 and n.GetTotalNumHs() == 1)
        if carbon.GetHybridization() == Chem.HybridizationType.SP3 and oh_count > 1:
            return False, "Extra hydroxyl groups found on core (might be phytosphingosine)"

    # 4. Check for long alkyl chains
    # One attached to the amide (fatty acyl) and one as the sphinganine tail
    alkyl_chain = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4]")
    chain_matches = len(mol.GetSubstructMatches(alkyl_chain))
    if chain_matches < 2:
        return False, "Missing required long alkyl chains"

    # 5. Basic composition checks
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 16:
        return False, f"Too few carbons ({c_count}) for N-acylsphinganine"
    if n_count != 1:
        return False, f"Must have exactly 1 nitrogen, found {n_count}"
    if o_count < 3:
        return False, f"Must have at least 3 oxygens, found {o_count}"

    # 6. Check for unsaturation in the core region
    # (to distinguish from sphingosine derivatives)
    double_bond = Chem.MolFromSmarts("C=C")
    if mol.HasSubstructMatch(double_bond):
        # Check if double bond is in the core region
        db_matches = mol.GetSubstructMatches(double_bond)
        core_atoms = set(mol.GetSubstructMatches(sphinganine_core)[0])
        for match in db_matches:
            if match[0] in core_atoms or match[1] in core_atoms:
                return False, "Contains double bond in core region (might be sphingosine derivative)"

    # The molecule passed all tests
    base_reason = "Valid N-acylsphinganine with correct core structure and stereochemistry"
    
    # Note if it's a glycosylated variant
    sugar_pattern = Chem.MolFromSmarts("[CH1,2]1[OH1,CH2][CH1][CH1][CH1][OH1,CH2]1")
    if mol.HasSubstructMatch(sugar_pattern):
        return True, base_reason + " (glycosylated variant)"
    
    return True, base_reason