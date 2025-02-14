"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: CHEBI:36326 lipopolysaccharide

Lipopolysaccharides (LPS) are complex natural compounds consisting of a trisaccharide repeating unit 
(two heptose units and octulosonic acid) with oligosaccharide side chains and 3-hydroxytetradecanoic 
acid units. They are a major constituent of the cell walls of Gram-negative bacteria.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.DataStructs import DiceSimilarity

def is_lipopolysaccharide(smiles):
    """
    Determines if a molecule is likely a lipopolysaccharide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a lipopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Check for presence of lipid chains
    lipid_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    has_lipid = mol.HasSubstructMatch(lipid_pattern)
    
    if not has_lipid:
        return False, "No lipid chains found"
    
    # Calculate topological fingerprint similarity to known LPS structures
    fps = [FingerprintMols.FingerprintMol(mol)]
    known_lps_mols = [Chem.MolFromSmiles(smi) for smi in known_lps_smiles]
    known_lps_fps = [FingerprintMols.FingerprintMol(m) for m in known_lps_mols]
    
    similarity_scores = [DiceSimilarity(fps[0], fp) for fp in known_lps_fps]
    avg_similarity = sum(similarity_scores) / len(similarity_scores)
    
    # LPS structures should have moderate similarity to known examples
    if avg_similarity < 0.4:
        return False, "Low structural similarity to known lipopolysaccharides"
    
    return True, "Molecule exhibits structural features characteristic of lipopolysaccharides"

# List of SMILES strings for known lipopolysaccharide structures
known_lps_smiles = [
    "O=C(OC(CCCC(=O)O[C@@H]1[C@@H](O)[C@H](O[C@@H]([C@H]1O)CO)O[C@H]2O[C@@H]([C@@H](O)[C@@H]([C@H]2O)O)CO)CCCCC)C(=O)O",
    "O=C(O[C@@H]1[C@@H](O[C@H](CO)[C@H]([C@@H]1OC(=O)CCCCCCCCCCCCCCC)OC(=O)C)OC[C@@H](O)[C@@H](O)CO)C(C)C",
    "CCCCCCCCCC(=O)O[C@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1O[C@@H](C[C@H](OC(C)=O)C(\C)=C\CO[C@H]([C@H](O)CO)[C@H](O)[C@H](O)CO)[C@](C)(O)CCC=C(C)C",
    "O=C(O[C@@H]1[C@@H](O[C@H](COC(=O)C)[C@H]([C@@H]1OC(=O)CCCCCCCCCCCCCCC)OC(=O)C)OC[C@@H](O)[C@@H](O)CO)CCCCC",
    "O=C(O[C@@H]1[C@@H](O[C@H](CO)[C@H]([C@@H]1OC(=O)CCCCCCCCCCCCC)OC(=O)C)OC[C@@H](O)[C@@H](O)CO)C(C)C"
]