import React, { ReactElement, useMemo } from "react";
import { TableTitle, TableTitleWrapper } from "../common/style";
import Table from "../Table";
import {
  aggregateCollectionsFromDatasets,
  useSourceData,
} from "src/common/queries/cellCards";

const Link = ({ title, url }: { title: string; url: string }) => {
  return (
    <a href={url} target="_blank">
      {title}
    </a>
  );
};
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

interface Collection {
  name: string;
  url: string;
  datasets: { id: string; label: string }[];
}

interface Collections {
  [name: string]: Collection;
}

const SourceDataTable = ({ cellTypeId }: Props) => {
  const { data } = useSourceData(cellTypeId);
  const { datasets = [] } = data;

  const collections: Collections = useMemo(() => {
    return aggregateCollectionsFromDatasets(datasets);
  }, [datasets]);

  const tableRows: TableRow[] = [];
  for (const collection of Object.keys(collections)) {
    tableRows.push({
      collection: (
        <Link
          title={collections[collection].name}
          url={collections[collection].url}
        />
      ),
      publication: "",
      tissue: "",
      disease: "",
      organism: "",
    });
  }

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
