"""
PDB Interface Analyzer

This module provides functionality to analyze protein-protein interfaces in PDB files.
Specifically, it identifies interface residues between two chains in a dimer structure.
"""

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
import numpy as np
from typing import Dict

def score_matrix(x,y) -> Dict[str, Dict[str, str]]:
    score_dict = {
        "A": {"A": "0", "R": "0", "N": "0", "D": "0", "C": "0", "Q": "0", "E": "0", "G": "0", "H": "0", "I": "+", "L": "+", "K": "0", "M": "0", "F": "0", "P": "0", "S": "0", "T": "0", "W": "0", "Y": "0", "V": "+"},
        "R": {"A": "0", "R": "--", "N": "0", "D": "++", "C": "0", "Q": "0", "E": "++", "G": "0", "H": "+", "I": "0", "L": "0", "K": "--", "M": "0", "F": "0", "P": "0", "S": "0", "T": "0", "W": "+", "Y": "+", "V": "0"},
        "N": {"A": "0", "R": "0", "N": "0", "D": "+", "C": "0", "Q": "+", "E": "0", "G": "0", "H": "+", "I": "0", "L": "0", "K": "0", "M": "0", "F": "0", "P": "0", "S": "+", "T": "+", "W": "0", "Y": "+", "V": "0"},
        "D": {"A": "0", "R": "++", "N": "+", "D": "--", "C": "0", "Q": "0", "E": "--", "G": "0", "H": "0", "I": "0", "L": "0", "K": "++", "M": "0", "F": "0", "P": "0", "S": "0", "T": "0", "W": "0", "Y": "0", "V": "0"},
        "C": {"A": "0", "R": "0", "N": "0", "D": "0", "C": "++", "Q": "0", "E": "0", "G": "0", "H": "0", "I": "0", "L": "0", "K": "0", "M": "0", "F": "0", "P": "0", "S": "0", "T": "0", "W": "0", "Y": "0", "V": "0"},
        "Q": {"A": "0", "R": "0", "N": "+", "D": "0", "C": "0", "Q": "0", "E": "+", "G": "0", "H": "+", "I": "0", "L": "0", "K": "0", "M": "0", "F": "0", "P": "0", "S": "0", "T": "0", "W": "0", "Y": "+", "V": "0"},
        "E": {"A": "0", "R": "++", "N": "0", "D": "--", "C": "0", "Q": "+", "E": "--", "G": "0", "H": "0", "I": "0", "L": "0", "K": "++", "M": "0", "F": "0", "P": "0", "S": "0", "T": "0", "W": "0", "Y": "0", "V": "0"},
        "G": {"A": "0", "R": "0", "N": "0", "D": "0", "C": "0", "Q": "0", "E": "0", "G": "0", "H": "0", "I": "0", "L": "0", "K": "0", "M": "0", "F": "0", "P": "0", "S": "0", "T": "0", "W": "0", "Y": "0", "V": "0"},
        "H": {"A": "0", "R": "+", "N": "+", "D": "0", "C": "0", "Q": "+", "E": "0", "G": "0", "H": "0", "I": "0", "L": "0", "K": "+", "M": "0", "F": "+", "P": "0", "S": "0", "T": "0", "W": "++", "Y": "++", "V": "0"},
        "I": {"A": "+", "R": "0", "N": "0", "D": "0", "C": "0", "Q": "0", "E": "0", "G": "0", "H": "0", "I": "++", "L": "++", "K": "0", "M": "+", "F": "+", "P": "-", "S": "0", "T": "0", "W": "+", "Y": "+", "V": "++"},
        "L": {"A": "+", "R": "0", "N": "0", "D": "0", "C": "0", "Q": "0", "E": "0", "G": "0", "H": "0", "I": "++", "L": "++", "K": "0", "M": "++", "F": "++", "P": "-", "S": "0", "T": "0", "W": "+", "Y": "+", "V": "++"},
        "K": {"A": "0", "R": "--", "N": "0", "D": "++", "C": "0", "Q": "0", "E": "++", "G": "0", "H": "+", "I": "0", "L": "0", "K": "--", "M": "0", "F": "0", "P": "0", "S": "0", "T": "0", "W": "0", "Y": "0", "V": "0"},
        "M": {"A": "0", "R": "0", "N": "0", "D": "0", "C": "0", "Q": "0", "E": "0", "G": "0", "H": "0", "I": "+", "L": "++", "K": "0", "M": "+", "F": "+", "P": "0", "S": "0", "T": "0", "W": "+", "Y": "+", "V": "+"},
        "F": {"A": "0", "R": "0", "N": "0", "D": "0", "C": "0", "Q": "0", "E": "0", "G": "0", "H": "+", "I": "+", "L": "++", "K": "0", "M": "+", "F": "++", "P": "-", "S": "0", "T": "0", "W": "++", "Y": "++", "V": "+"},
        "P": {"A": "0", "R": "0", "N": "0", "D": "0", "C": "0", "Q": "0", "E": "0", "G": "0", "H": "0", "I": "-", "L": "-", "K": "0", "M": "0", "F": "-", "P": "0", "S": "0", "T": "0", "W": "-", "Y": "-", "V": "-"},
        "S": {"A": "0", "R": "0", "N": "+", "D": "0", "C": "0", "Q": "0", "E": "0", "G": "0", "H": "0", "I": "0", "L": "0", "K": "0", "M": "0", "F": "0", "P": "0", "S": "0", "T": "0", "W": "0", "Y": "+", "V": "0"},
        "T": {"A": "0", "R": "0", "N": "+", "D": "0", "C": "0", "Q": "0", "E": "0", "G": "0", "H": "0", "I": "0", "L": "0", "K": "0", "M": "0", "F": "0", "P": "0", "S": "0", "T": "0", "W": "0", "Y": "+", "V": "0"},
        "W": {"A": "0", "R": "+", "N": "0", "D": "0", "C": "0", "Q": "0", "E": "0", "G": "0", "H": "++", "I": "+", "L": "+", "K": "0", "M": "+", "F": "++", "P": "-", "S": "0", "T": "0", "W": "++", "Y": "++", "V": "+"},
        "Y": {"A": "0", "R": "+", "N": "+", "D": "0", "C": "0", "Q": "+", "E": "0", "G": "0", "H": "++", "I": "+", "L": "+", "K": "0", "M": "+", "F": "++", "P": "-", "S": "+", "T": "+", "W": "++", "Y": "++", "V": "+"},
        "V": {"A": "+", "R": "0", "N": "0", "D": "0", "C": "0", "Q": "0", "E": "0", "G": "0", "H": "0", "I": "++", "L": "++", "K": "0", "M": "+", "F": "+", "P": "-", "S": "0", "T": "0", "W": "+", "Y": "+", "V": "++"}
    }
    score_mapping_dict = {"++": 2, "+": 1, "0": 0, "-": -1, "--": -2}
    score = score_mapping_dict[score_dict[x][y]]
    return score

def analyze_dimer_interface(pdb_file: str, distance_threshold: float = 5.0) -> Dict:
    """
    Analyze the interface between two chains in a PDB dimer structure.
    
    This function reads a PDB file containing exactly two amino acid chains (dimer),
    identifies the shorter chain, and finds all interface residues. An interface
    residue is defined as a residue in the shorter chain that has at least one
    atom within the distance threshold of any atom in the longer chain.
    
    Parameters
    ----------
    pdb_file : str
        Path to the PDB file containing the dimer structure
    distance_threshold : float, optional
        Maximum distance (in Angstroms) to consider residues as interacting
        Default is 5.0 Angstroms
    
    Returns
    -------
    dict
        A nested dictionary where:
        - Keys are tuples (chain_id, residue_number, residue_name) from the shorter chain
        - Values are dictionaries where:
            - Keys are tuples (chain_id, residue_number, residue_name) from the longer chain
            - Values are dictionaries containing {'distance': min_distance}
        
        Only residues from the shorter chain that have at least one contact
        within distance_threshold are included.
    
    Raises
    ------
    ValueError
        If the PDB file does not contain exactly 2 chains
        If the chains do not contain valid amino acids
    FileNotFoundError
        If the PDB file does not exist
    
    Examples
    --------
    >>> interface = analyze_dimer_interface('dimer.pdb')
    >>> # interface[('A', 25, 'LEU')][('B', 42, 'VAL')] = {'distance': 3.5}
    """
    # Parse PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('dimer', pdb_file)
    
    # Get all chains with amino acids
    chains = []
    for model in structure:
        for chain in model:
            # Check if chain has amino acid residues
            aa_residues = [res for res in chain if is_aa(res, standard=True)]
            if aa_residues:
                chains.append(chain)
    
    # Validate exactly 2 chains
    if len(chains) != 2:
        raise ValueError(f"PDB file must contain exactly 2 chains, found {len(chains)}")
    
    chain1, chain2 = chains[0], chains[1]
    
    # Get amino acid residues for each chain
    residues1 = [res for res in chain1 if is_aa(res, standard=True)]
    residues2 = [res for res in chain2 if is_aa(res, standard=True)]
    
    if not residues1 or not residues2:
        raise ValueError("Both chains must contain valid amino acid residues")
    
    # Determine which chain is shorter
    if len(residues1) <= len(residues2):
        short_chain_residues = residues1
        long_chain_residues = residues2
        short_chain = chain1
        long_chain = chain2
    else:
        short_chain_residues = residues2
        long_chain_residues = residues1
        short_chain = chain2
        long_chain = chain1
    
    # Build interface dictionary
    interface_dict = {}

    # Iterate through each residue in the shorter chain
    for short_res in short_chain_residues:
        short_key = (
            short_chain.id,
            short_res.id[1],  # residue number
            short_res.get_resname()
        )

        # Dictionary to store interacting residues from the long chain
        interacting_residues = {}

        # Check distance to all residues in the longer chain
        for long_res in long_chain_residues:
            # Calculate minimum distance between any atoms of the two residues
            min_distance = float('inf')

            for short_atom in short_res:
                for long_atom in long_res:
                    distance = np.linalg.norm(
                        short_atom.get_coord() - long_atom.get_coord()
                    )
                    if distance < min_distance:
                        min_distance = distance

            # If distance is within threshold, add to interacting residues
            if min_distance <= distance_threshold:
                long_key = (
                    long_chain.id,
                    long_res.id[1],  # residue number
                    long_res.get_resname()
                )
                # Get residue names (1-letter code) for scoring
                short_resname = short_res.get_resname()
                long_resname = long_res.get_resname()
                # Convert 3-letter to 1-letter code for score_matrix
                from Bio.Data.IUPACData import protein_letters_3to1
                try:
                    short_aa = protein_letters_3to1[short_resname.capitalize()]
                    long_aa = protein_letters_3to1[long_resname.capitalize()]
                    score = score_matrix(short_aa, long_aa)
                except Exception:
                    score = None
                interacting_residues[long_key] = {'distance': round(min_distance, 2), 'score': score}

        # Only add to interface dict if there are interacting residues
        if interacting_residues:
            interface_dict[short_key] = interacting_residues

    return interface_dict


def format_interface_report(interface_dict: Dict) -> str:
    """
    Format the interface dictionary into a human-readable report.
    
    Parameters
    ----------
    interface_dict : dict
        The interface dictionary returned by analyze_dimer_interface
    
    Returns
    -------
    str
        A formatted string report of the interface residues
    """
    if not interface_dict:
        return "No interface residues found."
    
    report_lines = [
        "=" * 70,
        "PDB DIMER INTERFACE ANALYSIS REPORT",
        "=" * 70,
        f"\nTotal interface residues in shorter chain: {len(interface_dict)}\n"
    ]
    
    for short_res_key, long_res_dict in sorted(interface_dict.items()):
        chain_id, res_num, res_name = short_res_key
        report_lines.append(f"\nResidue {res_name}-{res_num} (Chain {chain_id}):")
        report_lines.append(f"  Contacts with {len(long_res_dict)} residues in longer chain:")

        for long_res_key, dist_info in sorted(long_res_dict.items()):
            long_chain, long_num, long_name = long_res_key
            distance = dist_info['distance']
            score = dist_info.get('score', None)
            score_str = f", Score: {score}" if score is not None else ""
            report_lines.append(
                f"    - {long_name}-{long_num} (Chain {long_chain}): {distance:.2f} Å{score_str}"
            )
    
    report_lines.append("\n" + "=" * 70)
    return "\n".join(report_lines)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Analyze protein-protein interfaces in PDB dimer structures."
    )
    parser.add_argument(
        "pdb_file",
        type=str,
        help="Path to the PDB file containing the dimer structure"
    )
    parser.add_argument(
        "--distance",
        type=float,
        default=5.0,
        help="Distance threshold (in Angstroms) for defining interface residues (default: 5.0)"
    )
    
    args = parser.parse_args()
    
    try:
        interface = analyze_dimer_interface(args.pdb_file, args.distance)
        report = format_interface_report(interface)
        print(report)
    except Exception as e:
        print(f"Error: {e}")