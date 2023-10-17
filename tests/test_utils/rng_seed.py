import numpy as np


class SeededNumpyContext:
    def __init__(self, seed=0):
        self.seed = seed
        self.original_state = None

    def __enter__(self):
        self.original_state = np.random.get_state()
        np.random.seed(self.seed)

    def __exit__(self, exc_type, exc_val, exc_tb):
        np.random.set_state(self.original_state)
