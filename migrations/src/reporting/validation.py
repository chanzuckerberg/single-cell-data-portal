'''
Read validation errors from failing jobs. 
Validation errors result from running the cellxgene_schema.validate() method.
Validation errors are logged in the job-stream logs.
'''
from typing import List
import pandas as pd

_VALIDATION_COL_NAME = "validationError"

# filter log message for validation error messages
def _filter_logs_with_validation_message(log_results:List[str]):
    # only care about the messages with a valdiation message
    results = [lr for lr in log_results if lr.get("validation_message",None) is not None]

    # concatenate the error messages
    messages = [r["validation_message"] for r in results]
    return messages

# add a column to the dataframe to keep all validation errors
def merge_validation_errors(df:pd.DataFrame):
    if _VALIDATION_COL_NAME not in df.columns:
        df[_VALIDATION_COL_NAME] = df["logResults"].map(lambda lr: _filter_logs_with_validation_message(lr))