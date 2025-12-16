"""
Unit tests for the PDB Interface Analyzer
"""

import tempfile
import os
from pdb_interface_analyzer import analyze_dimer_interface


# Minimal PDB file content for testing (2 chains, simplified structure)
MINIMAL_DIMER_PDB = """ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00 20.00           C
ATOM      3  C   ALA A   1      11.500  11.000  10.000  1.00 20.00           C
ATOM      4  O   ALA A   1      11.000  12.000  10.000  1.00 20.00           O
ATOM      5  N   VAL A   2      12.000  10.000  11.000  1.00 20.00           N
ATOM      6  CA  VAL A   2      13.000  10.000  11.000  1.00 20.00           C
ATOM      7  C   VAL A   2      13.500  11.000  11.000  1.00 20.00           C
ATOM      8  O   VAL A   2      13.000  12.000  11.000  1.00 20.00           O
ATOM      9  N   LEU B   1      10.000  10.000  15.000  1.00 20.00           N
ATOM     10  CA  LEU B   1      11.000  10.000  15.000  1.00 20.00           C
ATOM     11  C   LEU B   1      11.500  11.000  15.000  1.00 20.00           C
ATOM     12  O   LEU B   1      11.000  12.000  15.000  1.00 20.00           O
ATOM     13  N   ILE B   2      12.000  10.000  16.000  1.00 20.00           N
ATOM     14  CA  ILE B   2      13.000  10.000  16.000  1.00 20.00           C
ATOM     15  C   ILE B   2      13.500  11.000  16.000  1.00 20.00           C
ATOM     16  O   ILE B   2      13.000  12.000  16.000  1.00 20.00           O
ATOM     17  N   GLY B   3      14.000  10.000  17.000  1.00 20.00           N
ATOM     18  CA  GLY B   3      15.000  10.000  17.000  1.00 20.00           C
ATOM     19  C   GLY B   3      15.500  11.000  17.000  1.00 20.00           C
ATOM     20  O   GLY B   3      15.000  12.000  17.000  1.00 20.00           O
END
"""

# PDB with only one chain (should raise error)
SINGLE_CHAIN_PDB = """ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00 20.00           C
ATOM      3  C   ALA A   1      11.500  11.000  10.000  1.00 20.00           C
ATOM      4  O   ALA A   1      11.000  12.000  10.000  1.00 20.00           O
END
"""

# PDB with three chains (should raise error)
THREE_CHAIN_PDB = """ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00 20.00           C
ATOM      3  N   VAL B   1      10.000  10.000  15.000  1.00 20.00           N
ATOM      4  CA  VAL B   1      11.000  10.000  15.000  1.00 20.00           C
ATOM      5  N   LEU C   1      10.000  10.000  20.000  1.00 20.00           N
ATOM      6  CA  LEU C   1      11.000  10.000  20.000  1.00 20.00           C
END
"""


def test_valid_dimer():
    """Test with a valid dimer structure"""
    print("Test 1: Valid dimer structure")
    print("-" * 70)
    
    # Create temporary PDB file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write(MINIMAL_DIMER_PDB)
        temp_pdb = f.name
    
    try:
        # Analyze the interface
        interface = analyze_dimer_interface(temp_pdb, distance_threshold=10.0)
        
        print(f"✓ Successfully analyzed dimer")
        print(f"  Found {len(interface)} interface residues")
        print(f"  Shorter chain has {len(interface)} residues at the interface")
        
        # Check that we got a dictionary
        assert isinstance(interface, dict), "Result should be a dictionary"
        
        # Check structure of returned data
        if interface:
            first_key = list(interface.keys())[0]
            assert isinstance(first_key, tuple), "Keys should be tuples"
            assert len(first_key) == 3, "Keys should have 3 elements (chain, resnum, resname)"
            
            first_value = interface[first_key]
            assert isinstance(first_value, dict), "Values should be dictionaries"
            
            if first_value:
                inner_key = list(first_value.keys())[0]
                assert isinstance(inner_key, tuple), "Inner keys should be tuples"
                assert len(inner_key) == 3, "Inner keys should have 3 elements"
                
                inner_value = first_value[inner_key]
                assert isinstance(inner_value, dict), "Inner values should be dictionaries"
                assert 'distance' in inner_value, "Inner dict should have 'distance' key"
                
                print(f"  Example: {first_key} → {inner_key}: {inner_value['distance']} Å")
        
        print("✓ Test passed!\n")
        return True
        
    except Exception as e:
        print(f"✗ Test failed: {e}\n")
        return False
    finally:
        os.unlink(temp_pdb)


def test_single_chain_error():
    """Test that single chain raises ValueError"""
    print("Test 2: Single chain (should raise error)")
    print("-" * 70)
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write(SINGLE_CHAIN_PDB)
        temp_pdb = f.name
    
    try:
        analyze_dimer_interface(temp_pdb)
        print("✗ Test failed: Should have raised ValueError\n")
        return False
    except ValueError as e:
        if "exactly 2 chains" in str(e):
            print(f"✓ Correctly raised ValueError: {e}")
            print("✓ Test passed!\n")
            return True
        else:
            print(f"✗ Test failed: Wrong error message: {e}\n")
            return False
    finally:
        os.unlink(temp_pdb)


def test_three_chain_error():
    """Test that three chains raises ValueError"""
    print("Test 3: Three chains (should raise error)")
    print("-" * 70)
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write(THREE_CHAIN_PDB)
        temp_pdb = f.name
    
    try:
        analyze_dimer_interface(temp_pdb)
        print("✗ Test failed: Should have raised ValueError\n")
        return False
    except ValueError as e:
        if "exactly 2 chains" in str(e):
            print(f"✓ Correctly raised ValueError: {e}")
            print("✓ Test passed!\n")
            return True
        else:
            print(f"✗ Test failed: Wrong error message: {e}\n")
            return False
    finally:
        os.unlink(temp_pdb)


def test_distance_threshold():
    """Test different distance thresholds"""
    print("Test 4: Distance threshold variations")
    print("-" * 70)
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write(MINIMAL_DIMER_PDB)
        temp_pdb = f.name
    
    try:
        # Test with large threshold (should find contacts)
        interface_large = analyze_dimer_interface(temp_pdb, distance_threshold=10.0)
        
        # Test with small threshold (likely fewer or no contacts)
        interface_small = analyze_dimer_interface(temp_pdb, distance_threshold=2.0)
        
        print(f"  With 10.0 Å threshold: {len(interface_large)} interface residues")
        print(f"  With 2.0 Å threshold: {len(interface_small)} interface residues")
        
        # Smaller threshold should have fewer or equal contacts
        assert len(interface_small) <= len(interface_large), \
            "Smaller threshold should not find more contacts"
        
        print("✓ Test passed!\n")
        return True
        
    except Exception as e:
        print(f"✗ Test failed: {e}\n")
        return False
    finally:
        os.unlink(temp_pdb)


def run_all_tests():
    """Run all tests"""
    print("=" * 70)
    print("RUNNING PDB INTERFACE ANALYZER TESTS")
    print("=" * 70)
    print()
    
    results = []
    results.append(("Valid dimer", test_valid_dimer()))
    results.append(("Single chain error", test_single_chain_error()))
    results.append(("Three chain error", test_three_chain_error()))
    results.append(("Distance threshold", test_distance_threshold()))
    
    print("=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✓ PASSED" if result else "✗ FAILED"
        print(f"{status}: {test_name}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    print("=" * 70)
    
    return passed == total


if __name__ == "__main__":
    success = run_all_tests()
    exit(0 if success else 1)
