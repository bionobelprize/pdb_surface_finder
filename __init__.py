"""
PDB Surface Finder - A tool for analyzing protein-protein interfaces

This package provides functionality to identify and analyze interface residues
in protein dimer structures from PDB files.
"""

from .pdb_interface_analyzer import analyze_dimer_interface, format_interface_report

__version__ = "1.0.0"
__all__ = ['analyze_dimer_interface', 'format_interface_report']
