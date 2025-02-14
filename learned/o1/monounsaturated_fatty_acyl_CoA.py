"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
from rdkit import Chem

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    A monounsaturated fatty acyl-CoA is any unsaturated fatty acyl-CoA in which the fatty acyl chain contains one carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify the thioester linkage: C(=O)S
    thioester_pattern = Chem.MolFromSmarts('C(=O)S')
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    
    if not thioester_matches:
        return False, "No thioester linkage found"
    
    # Process each thioester linkage
    for match in thioester_matches:
        carbonyl_c_idx = match[0]
        sulfur_idx = match[2]
        
        # Get the bond between carbonyl carbon and sulfur atom
        bond = mol.GetBondBetweenAtoms(carbonyl_c_idx, sulfur_idx)
        if bond is None:
            continue  # Bond not found, skip

        # Make an editable copy of the molecule
        mol_copy = Chem.RWMol(mol)
        # Remove the bond between carbonyl carbon and sulfur
        mol_copy.RemoveBond(carbonyl_c_idx, sulfur_idx)
        # Get the fragments resulting from bond removal
        frags = Chem.GetMolFrags(mol_copy.GetMol(), asMols=True, sanitizeFrags=True)
        
        # Identify the acyl chain fragment (contains the carbonyl carbon)
        acyl_chain_mol = None
        for frag in frags:
            atom_indices = [atom.GetIdx() for atom in frag.GetAtoms()]
            if carbonyl_c_idx in atom_indices:
                acyl_chain_mol = frag
                break

        if acyl_chain_mol is None:
            continue  # Acyl chain fragment not found, try next thioester

        # Analyze the acyl chain

        # Count the number of carbon atoms
        c_count = sum(1 for atom in acyl_chain_mol.GetAtoms() if atom.GetAtomicNum() == 6)
        # Count the number of heteroatoms (non-carbon, non-hydrogen atoms)
        heteroatom_count = sum(1 for atom in acyl_chain_mol.GetAtoms() if atom.GetAtomicNum() not in [1,6])

        # Check if acyl chain is composed mainly of carbons and hydrogens
        if heteroatom_count > 1:
            # Allow one oxygen atom (the carbonyl oxygen)
            return False, "Acyl chain contains heteroatoms, not a fatty acyl chain"

        # Check the length of the acyl chain (minimum 4 carbons for fatty acids)
        if c_count < 4:
            return False, f"Acyl chain too short ({c_count} carbons), not a fatty acyl chain"

        # Count the number of carbon-carbon double bonds
        double_bonds = 0
        for bond in acyl_chain_mol.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()
                # Count only carbon-carbon double bonds
                if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                    double_bonds +=1

        # Determine if the acyl chain is monounsaturated
        if double_bonds == 1:
            return True, "Contains one carbon-carbon double bond in the fatty acyl chain"
        elif double_bonds == 0:
            return False, "No carbon-carbon double bonds in the fatty acyl chain"
        else:
            return False, f"Contains {double_bonds} carbon-carbon double bonds in the fatty acyl chain"

    return False, "No suitable thioester linkage found"