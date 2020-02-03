from dcplib.config import Config


class LedgerDbConfig(Config):
    def __init__(self, *args, **kwargs):
        super().__init__(component_name='backend/ledger', secret_name='database', **kwargs)