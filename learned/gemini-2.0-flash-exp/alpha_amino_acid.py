"""
Classifies: CHEBI:33704 alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid has an amino group and a carboxylic acid group attached to the same carbon atom (alpha carbon).
    Handles zwitterions, substitutions and salts.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an alpha-amino acid, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove any counterions
    if "." in smiles:
        parts = smiles.split(".")
        if len(parts) > 1:
            # Try using the first part for core structure; this may not always work.
            mol = Chem.MolFromSmiles(parts[0])
            if mol is None:
              return False, "Invalid SMILES string after salt removal"
            
    # More flexible SMARTS pattern to capture variations in amino group and carboxylic acid.
    # [NX3,NX4+] allows both neutral and charged N. [CX3](=[OX1,OX2]) captures the carboxy carbon and variations in the oxygens.
    # [CX4] is the alpha carbon.
    # This pattern matches a carbon that is connected to both the amino and carboxyl group directly,
    # and it is used to define the alpha carbon relative to the other two
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4]([CX3](=[OX1,OX2])[OX1,OX2])")
    
    if not mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return False, "No alpha-amino acid substructure found (alpha C not connected)"
    
    # Count the number of carboxy groups and amino groups
    carboxy_pattern = Chem.MolFromSmarts("C(=O)[OX1,OX2]")
    amino_pattern = Chem.MolFromSmarts("[NX3,NX4+]")
    
    num_carboxy = len(mol.GetSubstructMatches(carboxy_pattern))
    num_amino = len(mol.GetSubstructMatches(amino_pattern))


    if num_carboxy != 1 or num_amino != 1:
       return False, f"Molecule does not contain exactly one carboxy group and one amino group: found {num_carboxy} carboxy groups and {num_amino} amino groups."
   
    # Molecular weight check (typical amino acid has MW < 300)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 300:
        return False, "Molecular weight is too high for a simple alpha-amino acid"

    return True, "Molecule matches the alpha-amino acid definition."