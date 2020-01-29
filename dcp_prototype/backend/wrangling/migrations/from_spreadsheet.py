from pandas import read_excel

from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.contributor import (
    Contributor,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.biosample_prep import (
    BiosamplePrep,
)

MAPPING_CLASSNAME_BY_SHEETNAME = {
    "Project - Contributors": Contributor,
    # "Project - Publications": Project,
    "Donor organism": BiosamplePrep,
}


class FromSpreadsheet:
    def __init__(self, spreadsheet_file_path):
        self.spreadsheet_file_path = spreadsheet_file_path

        spreadsheet_file = read_excel(spreadsheet_file_path, sheet_name=None)

        for sheet_name, sheet_contents in spreadsheet_file.items():
            print(sheet_name)
            if sheet_name in MAPPING_CLASSNAME_BY_SHEETNAME.keys():
                entity_class = MAPPING_CLASSNAME_BY_SHEETNAME[sheet_name]
                for index, row in sheet_contents.iterrows():
                    if index < 4:
                        continue
                    entity = entity_class.from_dcp_1_spreadsheet(row)
                    print(entity.name)
