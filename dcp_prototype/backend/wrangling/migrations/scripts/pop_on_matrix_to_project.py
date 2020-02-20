from argparse import ArgumentParser
import sys

sys.path.insert(0, "")  # noqa
import boto3
from botocore.exceptions import ClientError
from urllib.parse import urlparse
from dcp_prototype.backend.wrangling.migrations.utils.constants import SS2_QUANTIFICATION_PROTOCOL
from dcp_prototype.backend.wrangling.migrations.utils.id_generator import hca_accession_generator
from dcp_prototype.backend.ledger.code.common.ledger_orm import (
    DBSessionMaker,
    QuantificationProtocol,
    ExpressionFile,
    AnalysisFileQuantificationProtocolExpressionFileProcessJoin,
    AnalysisFile,
)
import os
import random

S3_CLIENT = boto3.client("s3")
S3_RESOURCE = boto3.resource("s3")


def upload_matrix_to_s3(matrix_file, s3_bucket):
    bucket = urlparse(s3_bucket).netloc
    prefix = urlparse(s3_bucket).path
    filename = matrix_file.split("/")[-1]

    file_prefix = prefix + filename
    full_s3_file_path = s3_bucket + filename

    try:
        S3_RESOURCE.Object(bucket, file_prefix).load()
        print(f"Matrix file s3://{file_prefix} already exists! Skipping upload.")
    except ClientError:
        S3_CLIENT.upload_file(matrix_file, bucket, file_prefix)
        print(f"Uploaded matrix file {file_prefix} to S3.")

    return full_s3_file_path, os.path.getsize(matrix_file)


def export_to_database(matrix_file_s3_uri, matrix_file_size):
    session_maker = DBSessionMaker()
    session = session_maker.session()

    # STEP 1: POPULATE ENTITIES

    # QuantificationProtocol table population
    quantification_protocols = {}
    quantification_protocol = SS2_QUANTIFICATION_PROTOCOL
    quantification_protocol.id = hca_accession_generator(QuantificationProtocol.__name__)
    quantification_protocols[quantification_protocol.id] = quantification_protocol
    session.add(quantification_protocol)

    # ExpressionFile table population
    expression_files = {}
    expression_file_id = hca_accession_generator(ExpressionFile.__name__)
    expression_file = ExpressionFile(
        id=expression_file_id,
        filename=matrix_file_s3_uri.split("/")[-1],
        file_format="loom",
        file_size=matrix_file_size,
        s3_uri=matrix_file_s3_uri,
    )
    expression_files[expression_file_id] = expression_file
    session.add(expression_file)

    # STEP 2: POPULATE JOIN TABLES

    # AnalysisFileQuantificationProtocolExpressionFileProcessJoin table population
    for analysis_file_entity in session.query(AnalysisFile):
        join_object_id = hca_accession_generator(AnalysisFileQuantificationProtocolExpressionFileProcessJoin.__name__)
        join_object = AnalysisFileQuantificationProtocolExpressionFileProcessJoin(
            id=join_object_id,
            analysis_file=analysis_file_entity,
            expression_file=random.choice(list(expression_files.values())),
            quantification_protocol=random.choice(list(quantification_protocols.values())),
        )
        session.add(join_object)

    print(f"Beginning commit to database")
    session.commit()
    print(f"Finished commit to database")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-i", "--input_matrix", nargs="+", required=False, help="A matrix to be uploaded to S3.",
    )
    parser.add_argument(
        "-b",
        "--bucket",
        nargs="+",
        required=False,
        help="The bucket into which files should be transfered if transferring local "
        "files as part of the project's metadata migration.",
    )

    arguments = parser.parse_args()

    if arguments.bucket:
        bucket = arguments.bucket[0]
        if "s3" in bucket:
            # Ensure that input directory ends in a slash
            if not bucket[-1] == "/":
                print(f"ERROR: S3 bucket must end in a slash!")
                sys.exit()

    input_file = arguments.input_matrix[0]
    bucket = arguments.bucket[0]

    # Upload matrix to bucket
    file_uri, file_size = upload_matrix_to_s3(input_file, bucket)

    # Upload new entities to Ledger DB
    export_to_database(file_uri, file_size)
