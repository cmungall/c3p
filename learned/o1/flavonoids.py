"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: flavonoids
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are characterized by a C6-C3-C6 skeleton consisting of two aromatic rings (A and B)
    connected via a heterocyclic ring (C) formed by the three-carbon bridge.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for different flavonoid subclasses

    # Flavone core: 2-phenylchromen-4-one
    flavone_pattern = Chem.MolFromSmarts('O=C1C=CC2=CC=CC=C2O1') 

    # Isoflavone core: 3-phenylchromen-4-one
    isoflavone_pattern = Chem.MolFromSmarts('O=C1C=CC(=CC1=O)c1ccccc1') 

    # Flavanone core: 2,3-dihydroflavone
    flavanone_pattern = Chem.MolFromSmarts('O=C1CC2=CC=CC=C2O1') 

    # Flavanol core: 3-hydroxyflavanone
    flavanol_pattern = Chem.MolFromSmarts('O[C@H]1COc2ccccc2C1=O') 

    # Chalcone core: 1,3-diphenylprop-2-en-1-one
    chalcone_pattern = Chem.MolFromSmarts('O=CC=CCc1ccccc1') 

    # Aurone core: 2-benzylidenebenzofuran-3(2H)-one
    aurone_pattern = Chem.MolFromSmarts('O=C1OC=CC2=CC=CC=C12') 

    # Anthocyanidin core
    anthocyanidin_pattern = Chem.MolFromSmarts('c1cc2c(c1)oc(=O)cc2c1ccccc1') 

    # General C6-C3-C6 skeleton with heterocyclic ring (ring C)
    flavonoid_core_pattern = Chem.MolFromSmarts('[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1-[#6]-[#6]-[#6]-2-[#6]~[#6]~[#6]~[#6]~[#6]~[#6]-2')

    # List of patterns with their corresponding flavonoid subclass names
    patterns = [
        ('Flavone', flavone_pattern),
        ('Isoflavone', isoflavone_pattern),
        ('Flavanone', flavanone_pattern),
        ('Flavanol', flavanol_pattern),
        ('Chalcone', chalcone_pattern),
        ('Aurone', aurone_pattern),
        ('Anthocyanidin', anthocyanidin_pattern),
        ('Flavonoid core', flavonoid_core_pattern),
    ]

    # Check for matches with flavonoid patterns
    for name, pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains {name} core structure"

    # Additional check for two aromatic rings connected via a heterocyclic ring
    ring_info = mol.GetRingInfo()
    aromatic_rings = [ring for ring in ring_info.AtomRings() if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]
    if len(aromatic_rings) >= 2:
        # Check for connection via a three-carbon bridge forming a heterocyclic ring
        hetero_ring_found = False
        for ring in ring_info.BondRings():
            atoms_in_ring = [mol.GetBondWithIdx(idx).GetBeginAtomIdx() for idx in ring] + \
                            [mol.GetBondWithIdx(idx).GetEndAtomIdx() for idx in ring]
            atoms_in_ring = set(atoms_in_ring)
            if len(atoms_in_ring) == 6 and any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in atoms_in_ring):
                hetero_ring_found = True
                break
        if hetero_ring_found:
            return True, "Contains two aromatic rings connected via a heterocyclic ring characteristic of flavonoids"

    return False, "No flavonoid core structure found"