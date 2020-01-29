from argparse import ArgumentParser
from from_spreadsheet import FromSpreadsheet

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_file",
        nargs=1,
        required=True,
        help="A XLSX file that was part of DCP 1.0 to be transformed into a valid DCP 2.0 spreadsheet.",
    )

    arguments = parser.parse_args()
    from_spreadsheet = FromSpreadsheet(arguments.input_file)


