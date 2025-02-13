"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: CHEBI:50254 nucleobase analogue
A molecule that can substitute for a normal nucleobase in nucleic acids.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for common nucleobase ring systems
    nucleobase_patterns = [
        Chem.MolFromSmarts("c1nc[nH]c1"),       # pyrimidine
        Chem.MolFromSmarts("c1ncnc1"),          # pyrazine
        Chem.MolFromSmarts("c1ncnc2[nH]cnc12"), # purine
        Chem.MolFromSmarts("c1ncnc2nc[nH]c12"), # imidazopyridine
    ]
    has_nucleobase_ring = any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns)

    if not has_nucleobase_ring:
        return False, "No nucleobase ring system found"

    # Check for common modifications
    modification_patterns = [
        Chem.MolFromSmarts("[NH2]"),            # amino
        Chem.MolFromSmarts("[OH]"),             # hydroxy
        Chem.MolFromSmarts("[OX2H]"),           # keto
        Chem.MolFromSmarts("[SX2]"),            # thiol
        Chem.MolFromSmarts("[NX3]"),            # azido
        Chem.MolFromSmarts("[NX2]=[OX1]"),      # nitro
        Chem.MolFromSmarts("[CH3]"),            # methyl
        Chem.MolFromSmarts("[CH2X4]"),          # alkyl
        Chem.MolFromSmarts("[F,Cl,Br,I]"),      # halogens
    ]
    has_modifications = any(mol.HasSubstructMatch(pattern) for pattern in modification_patterns)

    if not has_modifications:
        return False, "No modifications found on nucleobase ring"

    # Check for nucleic acid compatibility
    has_glycosidic_bond = mol.HasSubstructMatch(Chem.MolFromSmarts("[OX2][CX4]([CX4])[CX3]([OX2H])[CX3](=O)[NX3]"))
    if has_glycosidic_bond:
        return True, "Contains modified nucleobase ring and is compatible with nucleic acids"

    # Additional checks for specific functional groups
    has_primary_amine = mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]"))
    has_carbonyl = mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=[OX1])"))
    has_ring_nitrogen = any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms() if atom.IsInRing())

    if has_primary_amine and has_carbonyl and has_ring_nitrogen:
        return True, "Contains modified nucleobase ring with functional groups compatible for nucleic acid substitution"

    return False, "Modified nucleobase ring, but compatibility with nucleic acids is uncertain"