"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a non-classic icosanoid based on its SMILES string.
    A non-classic icosanoid is an oxygenated C20 fatty acid that is NOT a leukotriene or prostanoid.
    Lipoxins and resolvins are included.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple(bool, str): (True, reason) if molecule is a non-classic icosanoid, (False, reason) otherwise
                         or (None, None) if the SMILES string is invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
         return None, None

    # Check for 20 Carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count != 20:
        return False, f"Not a C20 fatty acid: {carbon_count} carbons"

    # Check for at least one oxygen
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count == 0:
        return False, "No oxygen atoms present"

     # Check for carboxylic acid
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group"

    # Check for long chain: 
    long_chain_pattern = Chem.MolFromSmarts("C[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Not a long-chain fatty acid"

    # Check for conjugated triene pattern (leukotriene-like)
    leukotriene_pattern_1 = Chem.MolFromSmarts("C=C-C=C-C=C")
    leukotriene_pattern_2 = Chem.MolFromSmarts("C=C-C=C-C=C-C")
    leukotriene_pattern_3 = Chem.MolFromSmarts("C=C-C=C-C=C-C=C")

    if mol.HasSubstructMatch(leukotriene_pattern_1) or \
       mol.HasSubstructMatch(leukotriene_pattern_2) or \
       mol.HasSubstructMatch(leukotriene_pattern_3):
       return False, "Likely a Leukotriene (conjugated triene pattern)"

    # Check for lipoxin/resolvin features - multiple double bonds and OH/epoxide/carbonyl
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    hydroxyl_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("[OH1]"))
    epoxide_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("C1OC1"))
    carbonyl_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("C=O"))

    if double_bond_count < 2 and len(hydroxyl_groups) == 1 and len(epoxide_groups) == 0 and len(carbonyl_groups) == 1:
       return False, "Too simple, likely a fatty acid."

    # Molecular Weight check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low"

    return True, "Likely a non-classic icosanoid"