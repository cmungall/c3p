"""
Classifies: CHEBI:29348 fatty amide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of an amide group (-C(=O)N-)
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"
    
    # Get all amide group matches
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) > 1: # check that only one amide group is present
        return False, "Too many amide groups present."

    # Check for a long carbon chain on the carbonyl side of the amide
    carbonyl_carbon_pattern = Chem.MolFromSmarts("[CX3](=[OX1])") #the carbonyl C
    carbonyl_carbon_matches = mol.GetSubstructMatches(carbonyl_carbon_pattern)

    #Check that amide carbonyls are attached to a long chain
    for match in carbonyl_carbon_matches:
        for atom_index in match:
            atom = mol.GetAtomWithIdx(atom_index)
            #find all neighbors of the carbonyl C
            neighbors = [neighbor.GetIdx() for neighbor in atom.GetNeighbors()]
            carbon_neighbors = [neighbor_idx for neighbor_idx in neighbors if mol.GetAtomWithIdx(neighbor_idx).GetAtomicNum() == 6]

            long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
            for neighbor_idx in carbon_neighbors:
                submol = Chem.PathToSubmol(mol, [atom_index, neighbor_idx])
                if submol.HasSubstructMatch(long_chain_pattern):
                    break
            else:
                  return False, "No long carbon chain attached to the amide carbonyl"


    # Check for additional carbonyl within the amide part
    additional_carbonyl = Chem.MolFromSmarts("N[CX3](=[OX1])[!#1]")
    if mol.HasSubstructMatch(additional_carbonyl):
        return False, "Additional carbonyl group found on the amide part."
      
    
    # Check molecular weight and number of carbons for long chain requirement
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 4: # At least 4 carbon chain
        return False, "Not enough carbons for fatty chain"
    if mol_wt < 100: # fatty amides should not be small, this is a very rough estimate.
        return False, "Molecular weight is too low for a typical fatty amide"

    return True, "Molecule is a fatty amide"