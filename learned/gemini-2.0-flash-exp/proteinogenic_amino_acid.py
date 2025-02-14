"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic alpha-amino acid structure, with allowance for labeled H atoms
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1:1][CX4:2]([CX3](=[OX1])O)([H,2H,3H])")
    match = mol.GetSubstructMatch(alpha_amino_acid_pattern)
    if not match:
          return False, "Not an alpha-amino acid structure"

    # Check L-configuration (except for glycine, which is achiral)
    alpha_carbon_idx = match[1]
    alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
    if alpha_carbon.GetDegree() == 4:  # chiral center, check the configuration
       stereo = alpha_carbon.GetChiralTag()

       if stereo == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
        return False, "Chirality is not L"
       if stereo == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
          pass #Correct L configuration
       else:
          return False, "No specified stereochemistry, must be L or achiral"


    # Check for the 20 standard amino acid sidechains (can be done via a list of SMARTS)
    # or allow any side chain
    # Check for selenocysteine (Sec)
    sec_pattern = Chem.MolFromSmarts("[NX3;H2,H1][CX4]([CX3](=[OX1])O)[CH2][SeH]")
    if mol.HasSubstructMatch(sec_pattern):
        return True, "Selenocysteine"
    
    # Check for pyrrolysine (Pyl)
    pyl_pattern = Chem.MolFromSmarts("C(=O)([C@@H](N)CCCCNC([C@H]1[C@@H](CC=N1)C)=O)O")
    if mol.HasSubstructMatch(pyl_pattern):
        return True, "Pyrrolysine"
    
    # Check for N-formylmethionine (fMet)
    fmet_pattern = Chem.MolFromSmarts("CNC(=O)[CH2][CH2][S][CH3]")
    if mol.HasSubstructMatch(fmet_pattern):
         return True, "N-Formylmethionine"
    
    # Passed all tests, the side chain is acceptable, and configuration is correct.
    return True, "Proteinogenic amino acid (side-chain check skipped)"