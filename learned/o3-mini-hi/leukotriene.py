"""
Classifies: CHEBI:25029 leukotriene
"""
"""
Classifies: Leukotrienes – any icosanoid stemming from a C20 polyunsaturated fatty acid derivative generated by leukocytes 
from arachidonic acid. They typically have a 20‐carbon acyclic chain with 4 double bonds (of which at least 3 are conjugated)
and a carboxyl (or related) moiety.
Improved heuristic:
  1) Molecule must be valid.
  2) It should have at most one ring.
  3) A carboxyl (or deprotonated carboxylate) group must be present.
  4) There must exist at least one contiguous (simple) subchain of exactly 20 carbon atoms that has 3 or 4 C=C bonds 
     (with at least 3 of them conjugated).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene derivative using improved heuristics.
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if molecule is classified as a leukotriene derivative.
       str: Explanation of the reasoning.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Rule 1: Check ring count: leukotrienes are largely acyclic (allow at most one small ring)
    ring_count = mol.GetRingInfo().NumRings()
    if ring_count > 1:
        return False, f"Too many rings ({ring_count}); leukotrienes are mostly acyclic (at most one small ring allowed)"
    
    # Rule 2: Check for a carboxyl group (or its deprotonated counterpart). 
    # This pattern covers carboxylic acids (C(=O)O) or carboxylates (C(=O)[O-]).
    carboxyl_smarts = Chem.MolFromSmarts("C(=O)[O;H1,-1]")
    if not mol.HasSubstructMatch(carboxyl_smarts):
        return False, "No carboxyl (or related) group found; expected in arachidonic acid derivatives"
    
    # Build carbon-only graph (nodes are carbon atom indices)
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_idxs) < 20:
        return False, f"Only {len(carbon_idxs)} carbon atoms present; a C20 subchain is required"
    
    # For each carbon atom, record its neighbors (only consider carbon-carbon bonds)
    carbon_neighbors = {idx: [] for idx in carbon_idxs}
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        # Include the edge if both atoms are carbons.
        if i in carbon_neighbors and j in carbon_neighbors:
            carbon_neighbors[i].append(j)
            carbon_neighbors[j].append(i)
    
    # We now perform a depth-first search (DFS) to look for any simple (non-repeating) path of exactly 20 carbon atoms.
    # For each candidate path, we check the bonds between consecutive carbons:
    #    Count the number of double bonds and the number of those that are conjugated.
    # Heuristic: Accept if there are either 3 or 4 double bonds, and at least 3 of them are conjugated.
    
    found_valid = False
    valid_explanation = ""
    
    # Use DFS with recursion. To avoid exponential blow-up in more complex systems,
    # we search only in molecules that typically aren’t huge.
    def dfs(path, visited):
        nonlocal found_valid, valid_explanation
        # If we already found a valid path, we can break out.
        if found_valid:
            return
        # When path length reaches 20, evaluate the bonds in the path.
        if len(path) == 20:
            dbl_count = 0
            conjugated_count = 0
            # Check bonds between consecutive carbons:
            for i in range(len(path) - 1):
                bond = mol.GetBondBetweenAtoms(path[i], path[i+1])
                if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                    dbl_count += 1
                    if bond.GetIsConjugated():
                        conjugated_count += 1
            # Accept only if double bond count is 3 or 4 and at least 3 are conjugated.
            if dbl_count in [3, 4] and conjugated_count >= 3:
                found_valid = True
                valid_explanation = f"Found a 20-carbon subchain with {dbl_count} C=C bonds ({conjugated_count} conjugated)"
            return
        
        # Continue DFS: iterate through neighbors of the last atom in the current path.
        current = path[-1]
        for nbr in carbon_neighbors.get(current, []):
            if nbr not in visited:
                visited.add(nbr)
                path.append(nbr)
                dfs(path, visited)
                if found_valid:
                    return
                path.pop()
                visited.remove(nbr)
    
    # Start DFS from each carbon atom.
    for start in carbon_idxs:
        dfs([start], set([start]))
        if found_valid:
            break
    
    if not found_valid:
        return False, "No appropriate 20-carbon subchain with the expected unsaturation pattern (3-4 C=C bonds with ≥3 conjugated) was found"
    
    # Report additional properties for explanation.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    explanation = (f"Molecule has a contiguous 20-carbon subchain where {valid_explanation}, "
                   f"a carboxyl group is present, MW is {mw:.1f} Da, and ring count is {ring_count} "
                   "– consistent with a leukotriene derivative.")
    return True, explanation
                        
# Example usage:
# Uncomment the following lines to test the function with a sample leukotriene A4 SMILES.
# example_smiles = "CCCCCC\\C=C/C\\C=C/C=C/C=C/[C@@H]1O[C@H]1CCCC(O)=O"
# result, reason = is_leukotriene(example_smiles)
# print(result, reason)