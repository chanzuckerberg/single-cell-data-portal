class PipelineStepMissing(ValueError):
    def __init__(self, array: str):
        super().__init__(f"Please run the pipeline step for generating {array} first.")
