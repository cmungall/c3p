"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is an aliphatic monocarboxylic acid with a chain length of less than C6 and
    no non-hydrocarbon substituents other than the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for non-hydrocarbon substituents (other than the carboxylic oxygens).
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num != 1 and atomic_num != 6 and atomic_num != 8:  # 1=H, 6=C, 8=O
            return False, "Non-hydrocarbon substituent found"


    # Count carbon atoms, including the carbonyl carbon, to check the chain length.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    if carbon_count > 6:
       return False, "More than 5 carbons in the chain"


    # Check for chain length (excluding carboxyl carbon).
    # For a molecule to be classified as a short chain, it must have a maximum of 5 carbons in the chain
    # However, this means the molecule can have 6 carbons, one at the carboxyl group and up to 5 in the aliphatic chain
    # Here we extract the carbon chain (excluding carbonyl carbon) and count its length.
    
    carbonyl_carbon_smarts = "[CX3](=[OX1])O"
    matches = mol.GetSubstructMatches(Chem.MolFromSmarts(carbonyl_carbon_smarts))
    if not matches:
      return False, "Could not locate carbonyl carbon, even though it was found before"

    
    carbon_chain_mol = Chem.RWMol(mol)
    carbon_chain_mol.RemoveAtom(matches[0][0]) #remove carbonyl carbon
    
    if len(carbon_chain_mol.GetAtoms()) > 5: #check the size of the new graph
        return False, "Too many carbons in the aliphatic chain"

    chain_c_count = sum(1 for atom in carbon_chain_mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if chain_c_count > 5:
        return False, "Too many carbons in the aliphatic chain"
    


    return True, "Short-chain fatty acid criteria met"