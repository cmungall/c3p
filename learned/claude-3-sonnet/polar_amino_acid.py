"""
Classifies: CHEBI:26167 polar amino acid
"""
"""
Classifies: CHEBI:27594 polar amino acid
Any amino acid whose side chain is capable of forming one or more hydrogen bonds.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondStereo

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polar amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for amino acid backbone pattern (N-C-C-C(=O)-O)
    backbone_pattern = Chem.MolFromSmarts("N[C@@H](C(=O)O)C")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No amino acid backbone found"
    
    # Look for common polar side chain groups
    polar_groups = ["[OH]",  # Hydroxyl
                    "[SH]",  # Thiol
                    "[NH2]",  # Primary amine
                    "[NH]",  # Secondary amine
                    "c[nH]c",  # Imidazole
                    "C(=O)N",  # Amide
                    "N=C(N)N"]  # Guanidino
    
    # Check if any polar group is present in the side chain
    for group in polar_groups:
        side_chain_pattern = Chem.MolFromSmarts(f"N[C@@H](C(=O)O)[C@@H]({group})")
        if mol.HasSubstructMatch(side_chain_pattern):
            return True, f"Contains {group} polar group in the side chain"
    
    # Check for zwitterionic forms (e.g., arginine)
    zwitterion_pattern = Chem.MolFromSmarts("N[C@@H](C(=O)[O-])[C@@H]([NH3+])")
    if mol.HasSubstructMatch(zwitterion_pattern):
        return True, "Zwitterionic form of a polar amino acid"
    
    return False, "No polar side chain group found"

# Examples
print(is_polar_amino_acid("N[C@@H](CC(O)=O)C(O)=O"))  # L-aspartic acid
print(is_polar_amino_acid("NCCCC[C@@H](N)C(O)=O"))  # L-lysine
print(is_polar_amino_acid("N[C@@H](CO)C(O)=O"))  # L-serine
print(is_polar_amino_acid("NC(CS)C(O)=O"))  # Cysteine
print(is_polar_amino_acid("N[C@@H](CC(N)=O)C(O)=O"))  # L-asparagine
print(is_polar_amino_acid("O=C(O)[C@H](N)CCC(=O)N"))  # D-glutamine
print(is_polar_amino_acid("CC(CCNC(N)=N)C(N)C(O)=O"))  # 3-methylarginine
print(is_polar_amino_acid("C[C@@H](O)[C@H](N)C(O)=O"))  # L-threonine
print(is_polar_amino_acid("N[C@@H](Cc1c[nH]c2ccccc12)C(O)=O"))  # L-tryptophan
print(is_polar_amino_acid("N[C@@H](CCC(N)=O)C(O)=O"))  # L-glutamine
print(is_polar_amino_acid("N[C@@H](Cc1c[nH]cn1)C(O)=O"))  # L-histidine
print(is_polar_amino_acid("NC(Cc1ccc(O)cc1)C(O)=O"))  # Tyrosine
print(is_polar_amino_acid("O=C([O-])[C@@H]([NH3+])[C@@H](CCNC(=[NH2+])N)C"))  # (3R)-3-methyl-L-arginine zwitterion