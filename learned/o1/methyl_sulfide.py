"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: methyl sulfide
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is an aliphatic sulfide where the sulfur atom is connected to two aliphatic carbon atoms,
    and at least one of those carbons is a methyl group. Peptides and larger biomolecules are excluded.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Exclude molecules with molecular weight above 700 Da (to exclude peptides)
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight > 700:
        return False, "Molecule too large (>{} Da) to be considered an aliphatic sulfide".format(int(mol_weight))

    # Identify sulfur atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Sulfur atom
            # Check if sulfur is connected to exactly two carbons
            neighbors = atom.GetNeighbors()
            if len(neighbors) != 2:
                continue
            carbon_neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 2:
                continue

            # Check that both carbons are aliphatic (sp3-hybridized and non-aromatic)
            aliphatic_carbons = []
            for carbon in carbon_neighbors:
                if carbon.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and not carbon.GetIsAromatic():
                    aliphatic_carbons.append(carbon)
            if len(aliphatic_carbons) != 2:
                continue

            # Check for peptide bond (amide linkage) involving these carbons (exclude peptides)
            is_peptide = False
            for carbon in aliphatic_carbons:
                for bond in carbon.GetBonds():
                    neighbor = bond.GetOtherAtom(carbon)
                    # Amide bond: carbonyl carbon (C=O) bonded to nitrogen (N)
                    if neighbor.GetAtomicNum() == 7:  # Nitrogen atom
                        for n_bond in neighbor.GetBonds():
                            n_neighbor = n_bond.GetOtherAtom(neighbor)
                            if n_neighbor.GetAtomicNum() == 6 and n_neighbor.GetIsAromatic() == False:
                                # Possible amide linkage
                                is_peptide = True
                                break
                if is_peptide:
                    break
            if is_peptide:
                continue

            # Check that at least one neighboring carbon is a methyl group
            methyl_found = False
            for carbon in aliphatic_carbons:
                if carbon.GetDegree() == 1:
                    # Carbon attached only to sulfur (methyl group)
                    methyl_found = True
                    break
            if not methyl_found:
                continue

            # Passed all checks
            return True, "Sulfur atom attached to a methyl group and another aliphatic carbon"
    
    return False, "Does not meet criteria for methyl sulfide (aliphatic sulfide with sulfur attached to methyl group)"