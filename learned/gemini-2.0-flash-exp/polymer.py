"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
import re

def is_polymer(smiles: str):
    """
    Determines if a molecule is likely a polymer based on its SMILES string.
    This is a heuristic approach and may not be 100% accurate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a polymer, False otherwise
        str: Reason for classification
    """

    if smiles is None or not isinstance(smiles,str) or len(smiles) == 0:
        return False, "Invalid SMILES string"

    # Handle mixtures and salts
    smiles_parts = smiles.split('.')
    if len(smiles_parts) > 1:
        is_polymer_mixture = False
        for part in smiles_parts:
             mol = Chem.MolFromSmiles(part)
             if mol is None:
                 return False, "Invalid SMILES string"
             mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
             num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
             if mol_wt > 500 and num_rotatable_bonds > 15:
                  chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]") # 8 carbon chain
                  if mol.HasSubstructMatch(chain_pattern):
                    is_polymer_mixture = True
                    break
             # Check for repeating isoprene units
             isoprene_pattern = Chem.MolFromSmarts("CC(C)=CC")
             if mol.HasSubstructMatch(isoprene_pattern) and mol_wt > 300 and num_rotatable_bonds > 10:
                 is_polymer_mixture = True
                 break
        if is_polymer_mixture:
             return True, "Mixture containing large molecule with repeating units or long chains, potentially a polymer"
        return False, "Mixture of small molecules, not a polymer"



    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"


    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

    # Check for long chains with high molecular weight and rotatable bonds
    if mol_wt > 1000 and num_rotatable_bonds > 30:
        return True, "High molecular weight and long chain, likely a polymer"
        
    # Check for repeating patterns
    chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]") # 8 carbon chain
    if mol.HasSubstructMatch(chain_pattern) and mol_wt > 500:
        return True, "Contains long carbon chain, likely a polymer"

    # Check for repeating isoprene units
    isoprene_pattern = Chem.MolFromSmarts("CC(C)=CC")
    if mol.HasSubstructMatch(isoprene_pattern) and mol_wt > 300 and num_rotatable_bonds > 10 :
            return True, "Contains repeating isoprene units, likely a polymer"

    # Check for phosphate, sulfate, carboxylate and other potential functional groups
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    sulfate_pattern = Chem.MolFromSmarts("S(=O)(=O)(O)")
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    ether_pattern = Chem.MolFromSmarts("C-O-C")
    ester_pattern = Chem.MolFromSmarts("C(=O)O")

    
    num_phosphates = len(mol.GetSubstructMatches(phosphate_pattern))
    num_sulfates = len(mol.GetSubstructMatches(sulfate_pattern))
    num_carboxylates = len(mol.GetSubstructMatches(carboxylate_pattern))
    num_amides = len(mol.GetSubstructMatches(amide_pattern))
    num_ethers = len(mol.GetSubstructMatches(ether_pattern))
    num_esters = len(mol.GetSubstructMatches(ester_pattern))

    if (num_phosphates + num_sulfates + num_carboxylates) > 2 and mol_wt > 500: #biopolymers
        return True, "Contains multiple phosphate, sulfate, or carboxylate groups, likely a biopolymer"
    if (num_amides > 4 or num_ethers > 5 or num_esters > 5) and mol_wt > 500:
      return True, "Contains multiple amides, ethers, or esters, potentially a polymer"

    # Check for explicit "poly" in the name.
    name = Chem.MolToInchiKey(mol)
    if name is not None and 'poly' in name.lower():
       return True, "Name contains poly prefix"
    
    #Check for oligomeric structures
    if mol_wt > 400 and num_rotatable_bonds > 10 and (num_amides > 2 or num_ethers > 2 or num_esters > 2):
        return True, "Contains repeating or long-chain structure, potentially an oligomer"

    return False, "Does not meet polymer criteria"