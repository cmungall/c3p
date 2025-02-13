"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: sulfolipid
A compound containing a sulfate ester group attached to a lipid-containing sugar moiety.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sulfate ester group (O-SO3H)
    sulfate_pattern = Chem.MolFromSmarts("[OX2][S](=[OX1])(=[OX1])[OX2H,OX1-]")
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "No sulfate ester group found"

    # Common sulfolipid sugar patterns
    galactosyl_pattern = Chem.MolFromSmarts("[CH2]1O[CH]([CH]([CH]([CH]([CH]1O)O)O)O)O")
    trehalose_pattern = Chem.MolFromSmarts("O[CH]1[CH](O)[CH](O)[CH]([CH](O1)CO)O")
    
    has_galactosyl = mol.HasSubstructMatch(galactosyl_pattern)
    has_trehalose = mol.HasSubstructMatch(trehalose_pattern)
    
    if not (has_galactosyl or has_trehalose):
        return False, "Missing characteristic sugar moiety (galactosyl or trehalose)"

    # Look for lipid characteristics
    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:  # Increased threshold for complex sulfolipids
        return False, "Too few carbons for a sulfolipid structure"

    # Look for long carbon chains (at least 12 carbons)
    long_chain = Chem.MolFromSmarts("[CH2,CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(long_chain):
        return False, "No long carbon chains found"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:  # Increased threshold for complex sulfolipids
        return False, "Molecular weight too low for sulfolipid"

    # Look for common structural features
    ester_pattern = Chem.MolFromSmarts("[#6]C(=O)O[#6]")
    amide_pattern = Chem.MolFromSmarts("[#6]C(=O)N[#6]")
    sphingosine_pattern = Chem.MolFromSmarts("[CH2]O[CH]([CH]([CH]([CH]O)O)O)[CH]([CH](O)/C=C/[#6])N")
    
    features = []
    if mol.HasSubstructMatch(ester_pattern):
        features.append("ester bonds")
    if mol.HasSubstructMatch(amide_pattern):
        features.append("amide bonds")
    if mol.HasSubstructMatch(sphingosine_pattern):
        features.append("sphingosine backbone")
    if has_galactosyl:
        features.append("galactosyl moiety")
    if has_trehalose:
        features.append("trehalose moiety")

    # Build reason string
    reason = "Contains sulfate ester group with characteristic sugar-lipid structure"
    if features:
        reason += f" including {', '.join(features)}"

    return True, reason