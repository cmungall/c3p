"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: Cannabinoids
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids are characterized by specific structural features such as
    the dibenzopyran ring system in phytocannabinoids (e.g., THC, CBD),
    long-chain polyunsaturated fatty acid derivatives linked to ethanolamine
    or glycerol (e.g., endocannabinoids like anandamide and 2-AG),
    and synthetic cannabinoids with varied core structures but common pharmacophores.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pattern for phytocannabinoids like THC, CBD (dibenzopyran core)
    dibenzopyran_pattern = Chem.MolFromSmarts('c1ccc2c(c1)Oc3ccccc3c2')  # General dibenzopyran ring

    # Check for phytocannabinoid features
    if mol.HasSubstructMatch(dibenzopyran_pattern):
        return True, "Contains dibenzopyran ring system characteristic of phytocannabinoids"

    # Define pattern for endocannabinoids like anandamide (N-acylethanolamines)
    n_acylethanolamine_pattern = Chem.MolFromSmarts('C(=O)NCCO')

    # Define pattern for 2-acylglycerols like 2-AG (monoacylglycerols)
    monoacylglycerol_pattern = Chem.MolFromSmarts('C(=O)OCC(O)CO')

    # Function to analyze fatty acid chains
    def _analyze_fatty_acid_chain(mol, carbonyl_carbon_idx):
        """
        Analyzes the fatty acid chain attached to the carbonyl carbon.

        Args:
            mol: RDKit Mol object
            carbonyl_carbon_idx: Index of the carbonyl carbon

        Returns:
            chain_length: Number of carbons in the chain
            num_double_bonds: Number of double bonds in the chain
        """
        atom = mol.GetAtomWithIdx(carbonyl_carbon_idx)
        neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        chain_start_idx = None
        for idx in neighbors:
            if idx != carbonyl_carbon_idx:
                chain_start_idx = idx
                break
        if chain_start_idx is None:
            return 0, 0

        visited = set()
        to_visit = [(chain_start_idx, None)]  # (atom_idx, bond)
        chain_length = 0
        num_double_bonds = 0

        while to_visit:
            current_idx, bond = to_visit.pop()
            if current_idx in visited:
                continue
            visited.add(current_idx)
            atom = mol.GetAtomWithIdx(current_idx)
            if atom.GetAtomicNum() != 6:
                continue
            chain_length += 1
            if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                num_double_bonds += 1
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:
                    b = mol.GetBondBetweenAtoms(current_idx, neighbor_idx)
                    to_visit.append((neighbor_idx, b))
        return chain_length, num_double_bonds

    # Check for endocannabinoid features (N-acylethanolamines)
    if mol.HasSubstructMatch(n_acylethanolamine_pattern):
        for match in mol.GetSubstructMatches(n_acylethanolamine_pattern):
            carbonyl_carbon_idx = match[0]  # Index of the carbonyl carbon
            chain_length, num_double_bonds = _analyze_fatty_acid_chain(mol, carbonyl_carbon_idx)
            if chain_length >= 16 and num_double_bonds >= 2:
                return True, "Contains long-chain polyunsaturated fatty acid amide linked to ethanolamine characteristic of endocannabinoids"
        return False, "Amide group found but fatty acid chain too short or not sufficiently unsaturated"

    # Check for endocannabinoid features (monoacylglycerols)
    if mol.HasSubstructMatch(monoacylglycerol_pattern):
        for match in mol.GetSubstructMatches(monoacylglycerol_pattern):
            carbonyl_carbon_idx = match[0]  # Index of the carbonyl carbon
            chain_length, num_double_bonds = _analyze_fatty_acid_chain(mol, carbonyl_carbon_idx)
            if chain_length >= 16 and num_double_bonds >= 2:
                return True, "Contains long-chain polyunsaturated fatty acid ester linked to glycerol characteristic of endocannabinoids"
        return False, "Ester group found but fatty acid chain too short or not sufficiently unsaturated"

    # Define pattern for synthetic cannabinoids (indole-based)
    indole_pattern = Chem.MolFromSmarts('c1ccc2c(c1)[nH]cc2')  # Indole core

    # Check for synthetic cannabinoid features (indole-based)
    if mol.HasSubstructMatch(indole_pattern):
        # Additional check for acyl or alkyl side chains
        side_chain_pattern = Chem.MolFromSmarts('n1cc(c2ccc(Cl)cc2)c(c1)C(=O)')
        if mol.HasSubstructMatch(side_chain_pattern):
            return True, "Contains indole core with acyl side chain characteristic of synthetic cannabinoids"
        return True, "Contains indole core characteristic of some synthetic cannabinoids"

    # Include pattern for naphthoylindole cannabinoids
    naphthoylindole_pattern = Chem.MolFromSmarts('c1ccc2c(c1)cccc2C(=O)N3CCCCC3')

    if mol.HasSubstructMatch(naphthoylindole_pattern):
        return True, "Contains naphthoylindole structure characteristic of synthetic cannabinoids"

    return False, "No characteristic cannabinoid structural features found"