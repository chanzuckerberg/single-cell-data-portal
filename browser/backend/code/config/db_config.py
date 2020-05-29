from dcplib.config import Config


class DataPortalDbConfig(Config):
    def __init__(self, *args, **kwargs):
        super().__init__(component_name="backend/data-portal", secret_name="database", **kwargs)
