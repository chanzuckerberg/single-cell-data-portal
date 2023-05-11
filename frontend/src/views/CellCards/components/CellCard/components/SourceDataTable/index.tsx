import React, { ReactElement, useMemo } from "react";
import { TableTitle, TableTitleWrapper } from "../common/style";
import Table from "../Table";
import Link from "../common/Link";
import { useSourceData } from "src/common/queries/cellCards";

interface TableRow {
  collection: ReactElement;
  publication: ReactElement;
  tissue: ReactElement;
  disease: ReactElement;
  organism: ReactElement;
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

const SourceDataTable = ({ cellTypeId }: Props) => {
  const { data: collections } = useSourceData(cellTypeId);

  const tableRows: TableRow[] = useMemo(() => {
    if (!collections) return [];
    const rows = [];
    for (const collection of collections) {
      const tissueNames = collection.tissue.map((tissue) => tissue.label);
      const diseaseNames = collection.disease.map((disease) => disease.label);
      const organismNames = collection.organism.map(
        (organism) => organism.label
      );

      rows.push({
        collection: (
          <Link
            title={collection.collection_name}
            url={collection.collection_url}
          />
        ),
        publication: (
          <Link
            title={collection.publication_title}
            url={`https://doi.org/${collection.publication_url}`}
          />
        ),
        tissue:
          tissueNames.length <= 7 ? (
            <div>
              {tissueNames.map((tissue) => {
                return <div>{tissue}</div>;
              })}
            </div>
          ) : (
            <div>{tissueNames.length} tissues</div>
          ),
        disease: (
          <div>
            {diseaseNames.map((disease) => {
              return <div>{disease}</div>;
            })}
          </div>
        ),
        organism: (
          <div>
            {organismNames.map((organism) => {
              return <div>{organism}</div>;
            })}
          </div>
        ),
      });
    }
    return rows;
  }, [collections]);

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
export default SourceDataTable;
