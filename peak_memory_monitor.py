#!/usr/bin/env python3
"""
Peak memory monitoring utilities for tracking memory usage during long operations.
"""

import logging
import threading
import time
from contextlib import contextmanager
from typing import Optional, Dict, Any

try:
    import psutil
    MEMORY_PROFILING_AVAILABLE = True
except ImportError:
    MEMORY_PROFILING_AVAILABLE = False


class PeakMemoryMonitor:
    """Monitor peak memory usage during operations."""
    
    def __init__(self, interval_seconds: float = 0.5):
        self.interval = interval_seconds
        self.monitoring = False
        self.peak_memory = 0
        self.peak_percent = 0
        self.samples = []
        self.monitor_thread = None
        self.process = psutil.Process() if MEMORY_PROFILING_AVAILABLE else None
        
    def start_monitoring(self, operation_name: str = "operation"):
        """Start background memory monitoring."""
        if not MEMORY_PROFILING_AVAILABLE:
            return
            
        self.monitoring = True
        self.peak_memory = 0
        self.peak_percent = 0
        self.samples = []
        self.operation_name = operation_name
        
        def monitor():
            while self.monitoring:
                try:
                    memory_info = self.process.memory_info()
                    memory_mb = memory_info.rss / 1024 / 1024
                    memory_percent = self.process.memory_percent()
                    
                    # Track peak
                    if memory_mb > self.peak_memory:
                        self.peak_memory = memory_mb
                        self.peak_percent = memory_percent
                    
                    # Store sample for analysis
                    self.samples.append({
                        'timestamp': time.time(),
                        'memory_mb': memory_mb,
                        'memory_percent': memory_percent
                    })
                    
                    # Keep only last 1000 samples to prevent memory buildup
                    if len(self.samples) > 1000:
                        self.samples = self.samples[-500:]
                        
                except Exception as e:
                    logging.warning(f"Memory monitoring error: {e}")
                
                time.sleep(self.interval)
        
        self.monitor_thread = threading.Thread(target=monitor, daemon=True)
        self.monitor_thread.start()
        logging.info(f"Started peak memory monitoring for {operation_name}")
    
    def stop_monitoring(self) -> Dict[str, Any]:
        """Stop monitoring and return peak memory stats."""
        if not MEMORY_PROFILING_AVAILABLE:
            return {}
            
        self.monitoring = False
        if self.monitor_thread:
            self.monitor_thread.join(timeout=2.0)
        
        stats = {
            'peak_memory_mb': self.peak_memory,
            'peak_memory_percent': self.peak_percent,
            'sample_count': len(self.samples),
            'duration_seconds': self.samples[-1]['timestamp'] - self.samples[0]['timestamp'] if self.samples else 0
        }
        
        logging.info(f"PEAK_MEMORY[{self.operation_name}]: {self.peak_memory:.1f} MB ({self.peak_percent:.1f}%) - {len(self.samples)} samples")
        
        return stats
    
    def get_memory_timeline(self) -> list:
        """Get the full memory timeline for plotting."""
        return self.samples.copy()


@contextmanager
def monitor_peak_memory(operation_name: str, interval_seconds: float = 0.5):
    """Context manager for monitoring peak memory during an operation."""
    monitor = PeakMemoryMonitor(interval_seconds)
    monitor.start_monitoring(operation_name)
    try:
        yield monitor
    finally:
        stats = monitor.stop_monitoring()


# Integration functions for existing codebase
def log_memory_usage_with_peak_monitoring(checkpoint_name: str = "", monitor: Optional[PeakMemoryMonitor] = None):
    """Enhanced memory logging that also shows current peak if monitoring."""
    if not MEMORY_PROFILING_AVAILABLE:
        return 0
        
    process = psutil.Process()
    memory_info = process.memory_info()
    memory_mb = memory_info.rss / 1024 / 1024
    memory_percent = process.memory_percent()
    
    if monitor and monitor.monitoring:
        peak_info = f" (Peak so far: {monitor.peak_memory:.1f} MB)"
    else:
        peak_info = ""
    
    logging.info(f"MEMORY_CHECKPOINT[{checkpoint_name}]: {memory_mb:.1f} MB RSS ({memory_percent:.1f}%){peak_info}")
    return memory_mb


# Example usage functions
def enhanced_tiledb_consolidate(array_path: str, ctx, **kwargs):
    """TileDB consolidate with peak memory monitoring."""
    import tiledb
    
    with monitor_peak_memory("TILEDB_CONSOLIDATE") as monitor:
        log_memory_usage_with_peak_monitoring("PRE_CONSOLIDATE", monitor)
        tiledb.consolidate(array_path, ctx=ctx, **kwargs)
        log_memory_usage_with_peak_monitoring("POST_CONSOLIDATE", monitor)


def enhanced_tiledb_vacuum(array_path: str, **kwargs):
    """TileDB vacuum with peak memory monitoring."""
    import tiledb
    
    with monitor_peak_memory("TILEDB_VACUUM") as monitor:
        log_memory_usage_with_peak_monitoring("PRE_VACUUM", monitor)
        tiledb.vacuum(array_path, **kwargs)
        log_memory_usage_with_peak_monitoring("POST_VACUUM", monitor)


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Simulate a memory-intensive operation
    with monitor_peak_memory("TEST_OPERATION") as monitor:
        log_memory_usage_with_peak_monitoring("START", monitor)
        
        # Simulate work
        import numpy as np
        data = np.random.random((10000, 1000))  # Allocate some memory
        time.sleep(2)
        
        log_memory_usage_with_peak_monitoring("MID", monitor)
        
        # Simulate more work
        data2 = np.random.random((20000, 1000))  # More memory
        time.sleep(2)
        
        log_memory_usage_with_peak_monitoring("END", monitor)
        
        del data, data2  # Cleanup
    
    print("Peak memory monitoring test completed!")