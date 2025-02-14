"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: CHEBI:37568 phosphatidic acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid is a derivative of glycerol where two hydroxy groups are esterified with fatty acids
    and the third is esterified with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define glycerol backbone pattern: three connected carbons each bonded to an oxygen
    glycerol_pattern = Chem.MolFromSmarts("[C;!R]([O])[C;!R]([O])[C;!R]([O])")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone with three oxygens found"

    # Define ester bond pattern: O-C(=O)
    ester_pattern = Chem.MolFromSmarts("[O][C](=O)[C]")
    # Define phosphate ester pattern: O-P(=O)(O)(O)
    phosphate_pattern = Chem.MolFromSmarts("[O][P](=O)(O)O")

    # Check each glycerol backbone match
    for match in glycerol_matches:
        c_atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        ester_count = 0
        phosphate_count = 0

        # Analyze oxygens attached to glycerol carbons
        for c_atom in c_atoms:
            o_atoms = [nbr for nbr in c_atom.GetNeighbors() if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in match]
            if not o_atoms:
                continue
            o_atom = o_atoms[0]  # Get the oxygen atom bonded to the glycerol carbon
            
            # Create a substructure starting from the oxygen
            sub_mol = Chem.FragmentOnBonds(mol, [mol.GetBondBetweenAtoms(c_atom.GetIdx(), o_atom.GetIdx()).GetIdx()])
            sub_mol = Chem.DeleteSubstructs(sub_mol, Chem.MolFromSmarts("[H]"))
            sub_smiles = Chem.MolToSmiles(sub_mol, isomericSmiles=True)

            # Check for ester linkage
            ester_mol = Chem.MolFromSmiles(sub_smiles)
            if ester_mol.HasSubstructMatch(ester_pattern):
                ester_count += 1
                continue

            # Check for phosphate linkage
            if ester_mol.HasSubstructMatch(phosphate_pattern):
                phosphate_count += 1
                continue

        if ester_count == 2 and phosphate_count == 1:
            return True, "Contains glycerol backbone with two fatty acid esters and one phosphoric acid ester"
    
    return False, "Does not match phosphatidic acid structure"