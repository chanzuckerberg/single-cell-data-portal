#!/usr/bin/env python3
"""
Memory profiling script for ATAC processing.

Usage:
    python memory_profile_atac.py <h5ad_file> <fragment_file>
    
    # With line-by-line profiling:
    python -m memory_profiler memory_profile_atac.py <h5ad_file> <fragment_file>
    
    # For detailed memory tracking:
    mprof run memory_profile_atac.py <h5ad_file> <fragment_file>
    mprof plot  # to visualize
"""

import logging
import sys
import tempfile
from pathlib import Path

# Set up logging to show memory checkpoints
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

def profile_atac_processing(h5ad_file: str, fragment_file: str):
    """Profile ATAC processing with the given files."""
    
    from backend.layers.processing.h5ad_data_file import H5ADDataFile
    from backend.layers.processing.utils.atac import log_memory_usage
    
    log_memory_usage("SCRIPT_START")
    
    print(f"Profiling ATAC processing for:")
    print(f"  H5AD file: {h5ad_file}")
    print(f"  Fragment file: {fragment_file}")
    
    try:
        # Load H5AD data
        log_memory_usage("PRE_H5AD_LOAD")
        h5ad_data_file = H5ADDataFile(h5ad_file)
        log_memory_usage("POST_H5AD_LOAD")
        
        # Create temporary CXG output directory
        with tempfile.TemporaryDirectory() as tmp_dir:
            cxg_output = Path(tmp_dir) / "test.cxg"
            
            log_memory_usage("PRE_CXG_CONVERSION")
            
            # Run the CXG conversion which includes ATAC processing
            h5ad_data_file.to_cxg(
                str(cxg_output), 
                sparse_threshold=25.0, 
                dataset_version_id="memory_profile_test",
                fragment_artifact_id=fragment_file
            )
            
            log_memory_usage("POST_CXG_CONVERSION")
            
    except Exception as e:
        print(f"Error during profiling: {e}")
        raise
    finally:
        log_memory_usage("SCRIPT_END")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python memory_profile_atac.py <h5ad_file> <fragment_file>")
        sys.exit(1)
        
    h5ad_file = sys.argv[1]
    fragment_file = sys.argv[2]
    
    # Verify files exist
    if not Path(h5ad_file).exists():
        print(f"H5AD file not found: {h5ad_file}")
        sys.exit(1)
        
    if not Path(fragment_file).exists():
        print(f"Fragment file not found: {fragment_file}")
        sys.exit(1)
    
    profile_atac_processing(h5ad_file, fragment_file)