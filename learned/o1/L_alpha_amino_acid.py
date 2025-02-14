"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
"""
Classifies: CHEBI:33709 L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.
    An L-alpha-amino acid has an alpha-carbon with L-configuration ('S' configuration,
    except for cysteine/selenocysteine which are 'R'), attached to an amino group and a carboxyl group.
    The molecule should not contain peptide bonds (amide linkages between amino acids).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens to the molecule to ensure accurate stereochemistry
    mol = Chem.AddHs(mol)

    # Assign stereochemistry to the molecule
    AllChem.AssignAtomChiralTagsFromStructure(mol)
    AllChem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # Define SMARTS patterns
    # Pattern for the alpha-amino acid backbone with chiral alpha-carbon
    alpha_amino_acid_pattern = Chem.MolFromSmarts('[NH2,NH3+][C@@H](*)C(=O)[O-,OH]')
    # Pattern for glycine (achiral)
    glycine_pattern = Chem.MolFromSmarts('NCC(=O)[O-,OH]')
    # Pattern to detect peptide bonds (exclude peptides)
    peptide_bond_pattern = Chem.MolFromSmarts('C(=O)N[C;!H0]')

    # Check for peptide bonds
    if mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "Molecule contains peptide bonds; not a single amino acid"

    # Check for alpha-amino acid backbone
    matches = mol.GetSubstructMatches(alpha_amino_acid_pattern)
    if not matches:
        # Check for glycine
        if mol.HasSubstructMatch(glycine_pattern):
            return True, "Molecule is glycine, an L-alpha-amino acid (achiral)"
        else:
            return False, "No alpha-amino acid backbone found"

    # Exclude molecules with multiple alpha-amino acid backbones
    if len(matches) > 1:
        return False, "Multiple alpha-amino acid units found; possible peptide"

    # Get the alpha-carbon atom
    match = matches[0]
    alpha_c_idx = match[1]
    alpha_c_atom = mol.GetAtomWithIdx(alpha_c_idx)

    # Check the CIP code of the alpha-carbon atom
    if alpha_c_atom.HasProp('_CIPCode'):
        cip_code = alpha_c_atom.GetProp('_CIPCode')

        # Identify side chain atoms (excluding amino and carboxyl groups)
        side_chain_atoms = [
            nbr.GetAtomicNum() for nbr in alpha_c_atom.GetNeighbors()
            if nbr.GetIdx() not in [match[0], match[2]] and nbr.GetAtomicNum() > 1
        ]

        has_sulfur = 16 in side_chain_atoms  # Sulfur atomic number
        has_selenium = 34 in side_chain_atoms  # Selenium atomic number

        if cip_code == 'S':
            if has_sulfur or has_selenium:
                return False, "Molecule is D-cysteine or D-selenocysteine (L-form has 'R' configuration)"
            else:
                return True, "Molecule is an L-alpha-amino acid with 'S' configuration at alpha-carbon"
        elif cip_code == 'R':
            if has_sulfur or has_selenium:
                return True, "Molecule is L-cysteine or L-selenocysteine with 'R' configuration at alpha-carbon"
            else:
                return False, "Molecule is D-alpha-amino acid"
        else:
            return False, "Unknown stereochemistry at alpha-carbon"
    else:
        # If alpha-carbon is not chiral, check if it's glycine
        neighbors = [
            nbr for nbr in alpha_c_atom.GetNeighbors()
            if nbr.GetIdx() not in [match[0], match[2]] and nbr.GetAtomicNum() > 1
        ]
        if len(neighbors) == 1 and neighbors[0].GetAtomicNum() == 1:
            # Only hydrogen as side chain, likely glycine
            return True, "Molecule is glycine, an L-alpha-amino acid (achiral)"
        else:
            return False, "Alpha-carbon has no stereochemistry assigned"

    return False, "No appropriate alpha-amino acid backbone with correct stereochemistry found"