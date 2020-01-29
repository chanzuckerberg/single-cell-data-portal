from pandas import read_excel

class FromSpreadsheet:
    def __init__(self, spreadsheet_file_path):
        self.spreadsheet_file_path = spreadsheet_file_path

        print(spreadsheet_file_path)
        spreadsheet_file = read_excel(open(spreadsheet_file_path, "rb"))
        print(spreadsheet_file.shape)