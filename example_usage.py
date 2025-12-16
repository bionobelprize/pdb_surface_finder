"""
Example usage of the PDB Interface Analyzer

This script demonstrates how to use the analyze_dimer_interface function
to analyze protein-protein interfaces in PDB files.
"""

from pdb_interface_analyzer import analyze_dimer_interface, format_interface_report


def main():
    """
    Example usage of the PDB interface analyzer.
    
    Note: You need to provide a valid PDB file path containing a dimer structure.
    """
    # Example 1: Basic usage with default parameters
    print("Example 1: Basic interface analysis")
    print("-" * 70)
    
    try:
        # Replace 'example_dimer.pdb' with your actual PDB file path
        pdb_file = 'example_dimer.pdb'
        
        # Analyze the dimer interface
        interface = analyze_dimer_interface(pdb_file)
        
        # Print formatted report
        report = format_interface_report(interface)
        print(report)
        
        # Access specific data
        print("\n\nExample of accessing specific data:")
        print("-" * 70)
        for short_res, contacts in list(interface.items())[:3]:  # Show first 3
            chain_id, res_num, res_name = short_res
            print(f"\n{res_name}-{res_num} (Chain {chain_id}) has {len(contacts)} contacts:")
            for long_res, dist_info in list(contacts.items())[:2]:  # Show first 2 contacts
                l_chain, l_num, l_name = long_res
                print(f"  → {l_name}-{l_num} (Chain {l_chain}): {dist_info['distance']} Å")
        
    except FileNotFoundError:
        print(f"Error: PDB file '{pdb_file}' not found.")
        print("Please provide a valid PDB file path.")
    except ValueError as e:
        print(f"Error: {e}")
    
    print("\n" + "=" * 70)
    
    # Example 2: Custom distance threshold
    print("\n\nExample 2: Custom distance threshold (3.5 Å)")
    print("-" * 70)
    
    try:
        pdb_file = 'example_dimer.pdb'
        
        # Use a stricter distance threshold
        interface = analyze_dimer_interface(pdb_file, distance_threshold=3.5)
        
        print(f"Interface residues with 3.5 Å cutoff: {len(interface)}")
        
    except (FileNotFoundError, ValueError) as e:
        print(f"Note: {e}")
    
    print("\n" + "=" * 70)
    
    # Example 3: Demonstrate data structure
    print("\n\nExample 3: Understanding the returned data structure")
    print("-" * 70)
    
    print("""
The function returns a nested dictionary with the following structure:

{
    (chain_id, residue_num, residue_name): {  # From shorter chain
        (chain_id, residue_num, residue_name): {  # From longer chain
            'distance': float  # Minimum distance in Angstroms
        },
        ...
    },
    ...
}

Example entry:
    interface[('A', 25, 'LEU')][('B', 42, 'VAL')] = {'distance': 3.5}

This means:
- LEU-25 in chain A (shorter chain) is an interface residue
- It contacts VAL-42 in chain B with a minimum distance of 3.5 Å
""")


if __name__ == "__main__":
    main()
