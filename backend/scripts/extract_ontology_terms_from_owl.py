from argparse import ArgumentParser
import boto3
import logging
from owlready2 import get_ontology


def extract_ontology_terms_from_file(input_file):
    """
    Returns the list of ontology classes in an OWL file.
    """

    logging.info(f"Reading in input OWL file {input_file}")
    ontology_file_object = open(input_file, "rb")
    ontology_object = get_ontology("")
    ontology_object.load(fileobj=ontology_file_object)

    ontology_classes = [ontology_class.name.replace("_", ":") for ontology_class in list(ontology_object.classes())]
    logging.info(f"Completed loading OWL file and found {len(ontology_classes)} ontology classes.")

    return ontology_classes


def generate_ontology_list_file(ontologies_list, filename, aws_bucket=None):
    """
    Generates a file containing all the ontology classes separated by a line break. If an AWS bucket is specified,
    uploaded the file to AWS S3.
    """

    file_contents = "\n".join(ontologies_list)

    if aws_bucket:
        logging.info(f"Writing ontology classes to S3 bucket {aws_bucket}.")

        # Truncate the filename to remove the path of the input file since the output filename is derived from the input
        # filename.
        truncated_filename = filename.split("/")[-1]

        s3_client = boto3.client("s3")
        s3_client.put_object(Body=file_contents, Bucket=aws_bucket, Key=truncated_filename)
        logging.info("Completed writing ontology classes to S3.")

    else:
        # Create a file locally
        logging.info(f"Writing ontology classes locally to {filename}.")

        output_file_object = open(filename, "w")
        output_file_object.write(file_contents)
        output_file_object.close()

        logging.info("Completed writing ontology classes to local file.")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_file",
        nargs="+",
        required=True,
        help="An OWL file that contains an ontology from which the classes will be extracted.",
    )
    parser.add_argument(
        "--upload_to_aws",
        action="store_true",
        help="If true, the resulting list of ontology classes extracted from the input OWL will be uploaded to AWS. An "
        "S3 bucket must be provided in order to do this via --aws_bucket.",
    )
    parser.add_argument(
        "-a",
        "--aws_bucket",
        nargs="+",
        help="The AWS bucket to which the ontology classes file that was generated should be uploaded.",
    )

    arguments = parser.parse_args()

    if arguments.input_file:
        input_file = arguments.input_file[0]

        # Verify that the input file is an OWL file
        if not input_file.endswith(".owl"):
            raise Exception(f"ERROR: Expecting an OWL file! Cannot parse {input_file}!")

    if arguments.upload_to_aws:
        # Verify that an AWS bucket has also been specified to which the ontologies file should be uploaded.
        if not arguments.aws_bucket:
            raise Exception(
                f"ERROR: Required to specify AWS bucket name if the ontologies file should be uploaded to "
                f"AWS via --aws_bucket."
            )

    if not arguments.upload_to_aws and arguments.aws_bucket:
        raise Exception(
            f"ERROR: Cannot specify an AWS bucket alone. Need to also add the --upload_to_aws flag to verify that you'd"
            f" like to upload to AWS."
        )

    if arguments.aws_bucket:
        aws_bucket = arguments.aws_bucket[0]

    list_of_ontologies = extract_ontology_terms_from_file(input_file)
    list_of_ontologies_filename = input_file.split(".")[0] + ".txt"

    if arguments.upload_to_aws:
        generate_ontology_list_file(list_of_ontologies, list_of_ontologies_filename, aws_bucket)
    else:
        generate_ontology_list_file(list_of_ontologies, list_of_ontologies_filename)
