"""
Classifies: CHEBI:35366 fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    Fatty acids are characterized by a long aliphatic chain and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
         return False, "Must have exactly one carboxylic acid group"


    # 2. Check for aliphatic nature (mostly C and H) and allowed substituents
    allowed_atoms = [1, 6, 7, 8, 9, 15, 16, 17, 35, 53]  # H, C, N, O, F, P, S, Cl, Br, I
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, "Non-allowed atoms present"
        
    # Count number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # 3. Check for chain length (min 4 carbons)
    if c_count < 4 :
        return False, "Chain too short to be a fatty acid"

    # 4. Branching Check - limit to 3 methyl or ethyl branches
    methyl_branch = Chem.MolFromSmarts("[CX4]([H])([H])([H])")
    ethyl_branch = Chem.MolFromSmarts("[CX4]([H])([H])([CX4][H])")
    branch_count = len(mol.GetSubstructMatches(methyl_branch)) + len(mol.GetSubstructMatches(ethyl_branch))
    if branch_count > 3:
        return False, "Too many branches for a typical fatty acid"

    # 5. Rings: check that if present, they are mostly made of carbons and that the carbon in the chain is not part of an "aromatic-like" cycle
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
         for ring in ring_info.AtomRings():
              is_embedded = True
              for atom_index in ring:
                    atom = mol.GetAtomWithIdx(atom_index)
                    if atom.GetAtomicNum() != 6:
                         is_embedded = False
                         break;
                    neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
                    carbon_count = 0
                    for n_idx in neighbors:
                        if mol.GetAtomWithIdx(n_idx).GetAtomicNum() == 6:
                             carbon_count += 1;
                    if carbon_count < 2: # if the carbon is not part of a chain, it's not an issue
                        is_embedded = False
                        break
              if is_embedded: # if it is a carbon ring
                    
                    is_aromatic = True
                    for atom_index in ring:
                          atom = mol.GetAtomWithIdx(atom_index)
                          if not atom.GetIsAromatic():
                               is_aromatic=False
                               break

                    if is_aromatic:
                        return False, "Aromatic ring found in the carbon chain"
    
    # 6. Check for unsaturation (double or triple bonds)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    has_double_bond = mol.HasSubstructMatch(double_bond_pattern)
    has_triple_bond = mol.HasSubstructMatch(triple_bond_pattern)

    # 7. Molecular weight check (minimum of 100, max of 1000)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
       return False, "Molecular weight too low to be a fatty acid"

    if mol_wt > 1000:
       return False, "Molecular weight too high to be a fatty acid"


    return True, "Meets criteria for a fatty acid"