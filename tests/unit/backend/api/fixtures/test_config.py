class TemporaryTestConfigChange:
    """
    Modify a backend.api.data_portal.config.app_config.*Config object while in the scope of this context manager, to
    support special test cases.
    """

    def __init__(self, config, config_props):
        self.orig_config_props = None
        self.tmp_config_props = config_props
        self.config = config

    def __enter__(self):
        self.orig_config_props = self.config.__dict__.copy()
        for k, v in self.tmp_config_props.items():
            self.config.__dict__[k] = v
        print(f"modified config for {self.config.__class__.__name__}")

    def __exit__(self, t, value, traceback):
        for k, v in self.orig_config_props.items():
            self.config.__dict__[k] = v
        print(f"restored config for {self.config.__class__.__name__}")
