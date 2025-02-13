"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: CHEBI:27772 mineral

A mineral is a chemical substance that is normally crystalline, formed as a result of geological processes.
It can include metamict substances, amorphous substances that have never been crystalline ('mineraloids'),
and biogenic minerals formed from geological processes on biogenic compounds.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for common mineral elements
    mineral_elements = [3, 4, 5, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 37, 38, 40, 42, 44, 47, 48, 49, 51, 52, 53, 55, 56, 57, 62, 63, 64, 65, 66, 67, 68, 70, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 90, 92]
    mol_elements = set(atom.GetAtomicNum() for atom in mol.GetAtoms())
    if not any(elem in mineral_elements for elem in mol_elements):
        return False, "Does not contain elements commonly found in minerals"
    
    # Check for common mineral ionic species
    ionic_species = ['[Cl-]', '[Br-]', '[I-]', '[F-]', '[OH-]', '[O-]', '[S-]', '[S--]', '[N-]', '[P-]', '[As-]', '[Li+]', '[Na+]', '[K+]', '[Rb+]', '[Cs+]', '[Be+2]', '[Mg+2]', '[Ca+2]', '[Sr+2]', '[Ba+2]', '[Al+3]', '[Fe+2]', '[Fe+3]', '[Zn+2]', '[Cu+]', '[Cu+2]', '[Ni+2]', '[Mn+2]', '[Cr+3]', '[Co+2]', '[NH4+]']
    mol_ions = Chem.MolToSmiles(mol, isomericSmiles=True).split('.')
    if any(ion in ionic_species for ion in mol_ions):
        return True, "Contains common mineral ionic species"
    
    # Check for common mineral bonding patterns
    ionic_pattern = Chem.MolFromSmarts("[+,+2,+3]~[-,-2]")
    covalent_pattern = Chem.MolFromSmarts("[O,S,Se,Te,N,P]~[O,S,Se,Te,N,P]")
    metallic_pattern = Chem.MolFromSmarts("[Fe,Cu,Ni,Co,Zn]=,@@[Fe,Cu,Ni,Co,Zn]")
    if mol.HasSubstructMatch(ionic_pattern) or mol.HasSubstructMatch(covalent_pattern) or mol.HasSubstructMatch(metallic_pattern):
        return True, "Contains common mineral bonding patterns"
    
    # Check for common mineral functional groups
    oxide_pattern = Chem.MolFromSmarts("[O-]=[#8,#16,#34,#52,#84]")
    halide_pattern = Chem.MolFromSmarts("[Cl,Br,I,F]~[Li,Na,K,Rb,Cs,Be,Mg,Ca,Sr,Ba,Al,Fe,Zn,Cu,Ni,Mn,Cr,Co]")
    sulfate_pattern = Chem.MolFromSmarts("[O-]S([O-])(=O)=O")
    phosphate_pattern = Chem.MolFromSmarts("[O-]P([O-])(=O)[O-]")
    carbonate_pattern = Chem.MolFromSmarts("[O-]C([O-])=O")
    if any(mol.HasSubstructMatch(pattern) for pattern in [oxide_pattern, halide_pattern, sulfate_pattern, phosphate_pattern, carbonate_pattern]):
        return True, "Contains common mineral functional groups"
    
    # Check for common mineral structural patterns
    layer_pattern = Chem.MolFromSmarts("[Si,Al]~[O]~[Si,Al]")
    chain_pattern = Chem.MolFromSmarts("[S,Se,Te]=[Fe,Cu,Ni,Co,Zn]~[S,Se,Te]=[Fe,Cu,Ni,Co,Zn]")
    if mol.HasSubstructMatch(layer_pattern) or mol.HasSubstructMatch(chain_pattern):
        return True, "Contains common mineral structural patterns"
    
    # If no mineral characteristics found, classify as non-mineral
    return False, "Does not exhibit typical mineral characteristics"