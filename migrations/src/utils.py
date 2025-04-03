from datetime import datetime,timezone

def dt_to_utc(dt:datetime) -> int:
    return int(dt.timestamp()*1_000)

def utc_to_dt(t:int) -> datetime:
    return datetime.fromtimestamp(int(t/1_000), timezone.utc) 

