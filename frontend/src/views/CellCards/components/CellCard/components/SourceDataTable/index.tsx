import React, { ReactElement } from "react";
import { TableTitle, TableTitleWrapper } from "../common/style";
import Table from "../Table";
import { useSourceData } from "src/common/queries/cellCards";

interface TableRow {
  collection: ReactElement;
  publication: string;
  tissue: string;
  disease: string;
  organism: string;
}
const tableColumns: Array<keyof TableRow> = [
  "collection",
  "publication",
  "tissue",
  "disease",
  "organism",
];

interface Props {
  cellTypeId: string;
}

const CanonicalMarkerGeneTable = ({ cellTypeId }: Props) => {
  const tableRows: TableRow[] = [];
  const { data, isLoading: _ } = useSourceData(cellTypeId);

  console.log(data);

  return (
    <>
      <TableTitleWrapper>
        <TableTitle>Data</TableTitle>
      </TableTitleWrapper>
      {tableRows.length && (
        <Table<TableRow> columns={tableColumns} rows={tableRows} />
      )}
    </>
  );
};
export default CanonicalMarkerGeneTable;
