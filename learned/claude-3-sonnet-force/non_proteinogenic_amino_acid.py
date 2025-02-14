"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
"""
Classifies: CHEBI:33709 non-proteinogenic amino acid

Any amino acid that is not naturally encoded in the genetic code of any organism.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a non-proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for core amino acid backbone (alpha-carbon with amino and carboxyl groups)
    amino_acid_pattern = Chem.MolFromSmarts("[C@H](N)(C(=O)O).*")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Missing core amino acid backbone"

    # Check for proteinogenic amino acids (exclude them)
    proteinogenic_smiles = ["N[C@@H](Cc1ccccc1)C(=O)O", "N[C@@H](Cc1c[nH]c2c1cccc2)C(=O)O",
                             "N[C@@H](CCc1ccccc1)C(=O)O", "N[C@@H](CO)C(=O)O",
                             "N[C@@H](CCS)C(=O)O", "N[C@@H](CCC(=O)N)C(=O)O",
                             "N[C@@H](CCC(=O)O)C(=O)O", "N[C@@H](Cc1c[nH]cn1)C(=O)O",
                             "N[C@@H](Cc1cnc[nH]1)C(=O)O", "N[C@@H](CC(C)C)C(=O)O",
                             "N[C@@H](CCC(N)=O)C(=O)O", "N[C@@H](CCCC(N)C(=O)O)C(=O)O",
                             "N[C@@H](CCSCC)C(=O)O", "N[C@@H](CC=C(C)C)C(=O)O",
                             "N[C@@H](CCC(N)=O)C(=O)O", "N[C@@H](CCc1c[nH]c2ccccc12)C(=O)O",
                             "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O", "N[C@@H](CS)C(=O)O",
                             "N[C@@H](CC1=CNC=N1)C(=O)O", "N[C@@H](CC1=CN=CN1)C(=O)O"]

    proteinogenic_mols = [Chem.MolFromSmiles(smi) for smi in proteinogenic_smiles]
    if any(mol.HasSubstructMatch(prot_mol) for prot_mol in proteinogenic_mols):
        return False, "Proteinogenic amino acid"

    # Check for additional modifications or functional groups
    has_mods = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [6, 7, 8, 9]:  # Allow only C, N, O, F
            has_mods = True
            break
        elif atom.GetIsAromatic():
            has_mods = True
            break
        elif atom.GetTotalDegree() > 4:
            has_mods = True
            break

    # Check molecular weight and atom counts
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    n_atoms = mol.GetNumAtoms()
    n_heavy_atoms = rdMolDescriptors.CalcNumHeavyAtoms(mol)

    if not has_mods and mol_wt < 100:
        return False, "Too small to be a non-proteinogenic amino acid"
    elif has_mods and mol_wt > 500:
        return False, "Molecular weight too high for an amino acid"
    elif n_atoms < 5 or n_heavy_atoms < 4:
        return False, "Too few atoms for an amino acid"

    return True, "Contains core amino acid backbone and additional modifications"