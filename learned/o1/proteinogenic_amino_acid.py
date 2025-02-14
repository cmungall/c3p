"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
"""
Classifies: CHEBI:59870 proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.
    A proteinogenic amino acid is one of the 23 alpha-amino acids that are incorporated into proteins during translation.
    Apart from glycine, which is non-chiral, all have L configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens to ensure proper valence
    mol = Chem.AddHs(mol)
    
    # Assign stereochemistry
    Chem.AssignAtomChiralTagsFromStructure(mol)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Define amino acid backbone SMARTS patterns
    # Corrected patterns with proper syntax
    backbone_smarts = '[N;H2][C@@H](*)C(=O)O'
    backbone = Chem.MolFromSmarts(backbone_smarts)
    if backbone is None:
        return False, "Invalid backbone SMARTS pattern"

    # Define glycine pattern (non-chiral)
    glycine_smarts = '[N;H2][CH2]C(=O)O' 
    glycine = Chem.MolFromSmarts(glycine_smarts)
    if glycine is None:
        return False, "Invalid glycine SMARTS pattern"

    # Check for glycine first
    if mol.HasSubstructMatch(glycine):
        return True, "Molecule is glycine, a proteinogenic amino acid"
    
    # Check for alpha-amino acid backbone
    matches = mol.GetSubstructMatches(backbone)
    if not matches:
        return False, "Molecule does not have alpha-amino acid backbone"
    
    # For each match, check chirality
    for match in matches:
        alpha_carbon_idx = match[1]  # Index of the alpha carbon
        alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
        chiral_tag = alpha_carbon.GetChiralTag()
        if chiral_tag == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            return False, "Alpha carbon is not chiral"
        
        # Get the stereochemistry (R or S)
        stereochemistry = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
        for idx, chiral_config in stereochemistry:
            if idx == alpha_carbon_idx:
                # For amino acids other than glycine, the alpha carbon should have L-configuration
                # Generally, L-amino acids have S configuration, except for cysteine and selenocysteine which are R
                side_chain_symbols = [nbr.GetSymbol() for nbr in alpha_carbon.GetNeighbors() 
                                      if nbr.GetAtomicNum() != 1 and nbr.GetIdx() not in match[:3]]
                if 'S' in side_chain_symbols or 'Se' in side_chain_symbols:
                    # For cysteine and selenocysteine, L-configuration corresponds to R
                    if chiral_config == 'R':
                        return True, "Molecule has L-configuration at alpha carbon, a proteinogenic amino acid"
                    else:
                        return False, "Cysteine and selenocysteine must have R configuration at alpha carbon"
                else:
                    # For all other amino acids, L-configuration corresponds to S
                    if chiral_config == 'S':
                        return True, "Molecule has L-configuration at alpha carbon, a proteinogenic amino acid"
                    else:
                        return False, "Amino acids must have L-configuration (S) at alpha carbon"
        return False, "Could not determine chirality"
    
    return False, "Molecule does not match proteinogenic amino acid criteria"