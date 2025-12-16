# PDB Surface Finder

A Python tool for analyzing protein-protein interfaces in PDB dimer structures.

## Overview

This tool reads PDB structure files containing protein dimers (two chains) and identifies interface residues - amino acids in the shorter chain that are within a specified distance threshold (default 5.0 Å) of amino acids in the longer chain.

## Features

- **Automatic Chain Identification**: Automatically determines which chain is shorter
- **Distance-based Interface Detection**: Finds residues within customizable distance thresholds
- **Detailed Interface Mapping**: Returns comprehensive data structure mapping all interface contacts
- **Input Validation**: Ensures the PDB file contains exactly 2 chains
- **Flexible Output**: Provides both raw data and formatted reports

## Installation

### Requirements

- Python 3.7+
- BioPython
- NumPy

### Install Dependencies

```bash
pip install -r requirements.txt
```

## Usage

### Basic Usage

```python
from pdb_interface_analyzer import analyze_dimer_interface, format_interface_report

# Analyze a dimer PDB file
interface = analyze_dimer_interface('your_dimer.pdb')

# Print formatted report
report = format_interface_report(interface)
print(report)
```

### Custom Distance Threshold

```python
# Use a stricter 3.5 Å cutoff
interface = analyze_dimer_interface('your_dimer.pdb', distance_threshold=3.5)
```

### Accessing Interface Data

The function returns a nested dictionary structure:

```python
interface = analyze_dimer_interface('dimer.pdb')

# Structure: {short_chain_residue: {long_chain_residue: {'distance': float}}}
# Example access:
for short_res, contacts in interface.items():
    chain_id, res_num, res_name = short_res
    print(f"Residue {res_name}-{res_num} (Chain {chain_id}):")
    
    for long_res, dist_info in contacts.items():
        l_chain, l_num, l_name = long_res
        distance = dist_info['distance']
        print(f"  → {l_name}-{l_num} (Chain {l_chain}): {distance} Å")
```

## Data Structure

The returned dictionary has the following structure:

```
{
    (chain_id, residue_num, residue_name): {  # From shorter chain
        (chain_id, residue_num, residue_name): {  # From longer chain
            'distance': float  # Minimum distance in Angstroms
        },
        ...
    },
    ...
}
```

**Example:**
```python
interface[('A', 25, 'LEU')][('B', 42, 'VAL')] = {'distance': 3.5}
```

This means:
- LEU-25 in chain A (shorter chain) is an interface residue
- It contacts VAL-42 in chain B with a minimum distance of 3.5 Å

## Function Reference

### `analyze_dimer_interface(pdb_file, distance_threshold=5.0)`

Analyze the interface between two chains in a PDB dimer structure.

**Parameters:**
- `pdb_file` (str): Path to the PDB file containing the dimer structure
- `distance_threshold` (float, optional): Maximum distance (in Angstroms) to consider residues as interacting. Default is 5.0 Å

**Returns:**
- `dict`: Nested dictionary of interface residues and their contacts

**Raises:**
- `ValueError`: If the PDB file does not contain exactly 2 chains
- `FileNotFoundError`: If the PDB file does not exist

### `format_interface_report(interface_dict)`

Format the interface dictionary into a human-readable report.

**Parameters:**
- `interface_dict` (dict): The interface dictionary returned by analyze_dimer_interface

**Returns:**
- `str`: A formatted string report of the interface residues

## Testing

Run the test suite:

```bash
python test_pdb_interface.py
```

## Examples

See `example_usage.py` for detailed examples of how to use the analyzer.

## Algorithm

1. Parse the PDB file and extract all chains containing amino acids
2. Validate that exactly 2 chains are present
3. Identify the shorter chain
4. For each residue in the shorter chain:
   - Calculate the minimum distance to all atoms in each residue of the longer chain
   - If the minimum distance is ≤ threshold, mark it as an interface residue
   - Store the contact information with distance
5. Return the complete interface dictionary

## License

This project is open source.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.