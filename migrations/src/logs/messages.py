'''
Utility functions for parsing logger information stored with each job
'''
from typing import List
import json

# get all log messages at ERROR log level 
def filter_errors(events):
    return [e for e in events if e["levelname"]=="ERROR"]

# get all log messages at INFO log level 
def filter_info(events):
    return [e for e in events if e["levelname"]=="INFO"]

# parse error message info from a specific message
def extract_error(event):
    keys = ["message","lineno","pathname","exc_info"]
    return {k:event.get(k,None) for k in keys}

