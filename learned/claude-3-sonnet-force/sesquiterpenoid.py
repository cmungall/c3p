"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies: CHEBI:36198 Sesquiterpenoid

A sesquiterpenoid is any terpenoid derived from a sesquiterpene.
The term includes compounds in which the C15 skeleton of the parent sesquiterpene 
has been rearranged or modified by the removal of one or more skeletal atoms 
(generally methyl groups).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for C15 skeleton
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 15:
        return False, "Not a C15 skeleton"
    
    # Check for terpene-like structure
    ring_info = mol.GetRingInfo()
    if not any(len(ring) >= 6 for ring in ring_info.AtomRings()):
        return False, "No cyclohexane-like ring found (typical of terpenes)"
    
    # Check for common sesquiterpene scaffolds
    scaffolds = ["C[C@@H]1CCC2=CC(=O)CC[C@]12C", # germacrane
                 "CC1=CCC2C(C1)C1=C(CC[C@@]3(C)C[C@@H](OC3=O)C1)C2", # eudesmane
                 "CC1=C2C(C=C(C)C)C1=CCC1C2CCC2=CC(=O)CC12", # guaiane
                 "CC1=CCC2C(C1)C1=C3CC[C@@H](C)C(C1)C3CCC2=O", # cadinane
                 "CC1=CCC2C(CCC3CC2=CC(=O)C23C)C1", # muurolane
                 "CC1=CCC2C3(C)CCC(C)(C3CCC2C1=O)C", # longibornane
                 "CC1=C2CCC3C(=CCC4CC(CC(C4)C3=C)C2=C)C1", # himachalane
                 "CC(C)C1=CCC=C2C=CC(C)(CCC3CC2=CC(C31)C)C", # spiroaxane
                 "CC1CCC2=C3CCC4=CC(=O)C=CC4=C3CCC12C", # africane
                 "CC1=C2CC[C@@H](C)C(=O)C2=CC2=C1CC[C@@H](O)C2" # cedrane
                ]
    
    for scaffold in scaffolds:
        scaffold_mol = Chem.MolFromSmarts(scaffold)
        if mol.HasSubstructMatch(scaffold_mol):
            return True, f"Contains {scaffold_mol.GetProp('_Name')} sesquiterpene scaffold"
    
    # Check for structural rearrangements
    rearrangements = ["[C@@H]1CC=C(C)C(=O)C1", # cyclobutanone
                      "[CH2]C(=O)C", # ketone
                      "C=CC(=O)O", # carboxylic acid
                      "OC" # ether/alcohol
                     ]
    mod_count = 0
    for rearrange in rearrangements:
        rearrange_mol = Chem.MolFromSmarts(rearrange)
        mod_count += len(mol.GetSubstructMatches(rearrange_mol))
    
    if mod_count > 2:
        return True, "C15 scaffold modified by rearrangements/functional groups"
    
    return False, "Could not confidently classify as a sesquiterpenoid"