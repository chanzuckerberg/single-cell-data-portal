#!/usr/bin/env python3
"""
Advanced memory tracing for ATAC processing using tracemalloc.
"""

import tracemalloc
import logging
from typing import Optional

def detailed_memory_trace(h5ad_file: str, fragment_file: str):
    """Run ATAC processing with detailed memory tracing."""
    
    # Start tracing
    tracemalloc.start()
    
    # Take initial snapshot
    snapshot1 = tracemalloc.take_snapshot()
    
    try:
        from backend.layers.processing.h5ad_data_file import H5ADDataFile
        
        print(f"Starting ATAC processing with tracemalloc...")
        
        h5ad_data_file = H5ADDataFile(h5ad_file)
        
        # Take snapshot after H5AD load
        snapshot2 = tracemalloc.take_snapshot()
        top_stats = snapshot2.compare_to(snapshot1, 'lineno')
        
        print("\n=== Memory usage after H5AD load ===")
        for stat in top_stats[:5]:
            print(stat)
        
        # Continue with CXG conversion
        with tempfile.TemporaryDirectory() as tmp_dir:
            cxg_output = Path(tmp_dir) / "test.cxg"
            
            h5ad_data_file.to_cxg(
                str(cxg_output),
                sparse_threshold=25.0,
                dataset_version_id="trace_test", 
                fragment_artifact_id=fragment_file
            )
        
        # Final snapshot
        snapshot3 = tracemalloc.take_snapshot()
        top_stats = snapshot3.compare_to(snapshot2, 'lineno')
        
        print("\n=== Memory usage after CXG conversion ===")
        for stat in top_stats[:10]:
            print(stat)
            
        print(f"\n=== Peak memory usage ===")
        current, peak = tracemalloc.get_traced_memory()
        print(f"Current: {current / 1024 / 1024:.1f} MB")
        print(f"Peak: {peak / 1024 / 1024:.1f} MB")
        
    finally:
        tracemalloc.stop()


if __name__ == "__main__":
    import sys
    import tempfile
    from pathlib import Path
    
    if len(sys.argv) != 3:
        print("Usage: python memory_trace_atac.py <h5ad_file> <fragment_file>")
        sys.exit(1)
        
    detailed_memory_trace(sys.argv[1], sys.argv[2])