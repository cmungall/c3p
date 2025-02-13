"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: tetrachlorobenzene
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule contains a tetrachlorobenzene substructure.
    Tetrachlorobenzene is a benzene ring with exactly 4 chlorine atoms attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all benzene rings
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    
    if not benzene_matches:
        return False, "No benzene ring found"
    
    # For each benzene ring, check if it has exactly 4 chlorines attached
    for benzene_match in benzene_matches:
        benzene_atoms = set(benzene_match)
        chlorine_count = 0
        
        # Check each carbon in the benzene ring for attached chlorines
        for atom_idx in benzene_match:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 17:  # Chlorine atomic number
                    chlorine_count += 1
        
        if chlorine_count == 4:
            return True, "Found benzene ring with exactly 4 chlorine atoms"
            
    return False, "No benzene ring with exactly 4 chlorines found"

def test_examples():
    """Test function with known examples"""
    examples = [
        "Clc1cc(Cl)c(Cl)cc1Cl",  # 1,2,4,5-tetrachlorobenzene
        "Clc1ccc(Cl)c(Cl)c1Cl",  # 1,2,3,4-tetrachlorobenzene
        "Clc1cc(Cl)c(Cl)c(Cl)c1",  # 1,2,3,5-tetrachlorobenzene
        "CC1=CC=CC=C1",  # Toluene (negative example)
        "Clc1c(Cl)c(Cl)c(Cl)c(Cl)c1Cl"  # Hexachlorobenzene (negative example)
    ]
    
    for smiles in examples:
        result, reason = is_tetrachlorobenzene(smiles)
        print(f"SMILES: {smiles}")
        print(f"Result: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_examples()