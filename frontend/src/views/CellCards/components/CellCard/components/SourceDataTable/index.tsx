import React, { ReactElement, useState, useEffect } from "react";
import { SOURCE_DATA_SECTION_ID } from "../../../CellCardSidebar";
import { TableTitle, TableTitleWrapper } from "../common/style";
import Table from "../Table";

// import {
//   aggregateCollectionsFromDatasets,
//   useSourceData,
// } from "src/common/queries/cellCards";

const Link = ({ title, url }: { title: string; url: string }) => {
  return (
    <a href={url} target="_blank">
      {title}
    </a>
  );
};
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

interface Collection {
  collection_name: string;
  collection_url: string;
  publication_url: string;
  publication_title: string;
  tissue: { label: string; ontology_term_id: string }[];
  disease: { label: string; ontology_term_id: string }[];
  organism: { label: string; ontology_term_id: string }[];
}

interface Props {
  cellTypeId: string;
}

// interface Collection {
//   name: string;
//   url: string;
//   datasets: { id: string; label: string }[];
// }

// interface Collections {
//   [name: string]: Collection;
// }

const MAX_RETRIES = 3; // Maximum number of retries

async function fetchData(
  cellTypeId: string,
  retries: number = MAX_RETRIES
): Promise<any> {
  try {
    const res = await fetch(`/api/sourceData?cellTypeId=${cellTypeId}`);

    if (!res.ok) {
      throw new Error(res.statusText);
    }

    const data = await res.json();
    return data;
  } catch (error) {
    if (retries > 0) {
      return fetchData(cellTypeId, retries - 1);
    } else {
      throw new Error("Fetching data failed after several retries");
    }
  }
}

const SourceDataTable = ({ cellTypeId }: Props) => {
  const [collections, setCollections] = useState<Collection[]>([]);
  // const { data } = useSourceData(cellTypeId);
  // const { datasets = [] } = data;

  // const collections: Collections = useMemo(() => {
  //   return aggregateCollectionsFromDatasets(datasets);
  // }, [datasets]);
  useEffect(() => {
    fetchData(cellTypeId)
      .then((data) => setCollections(data))
      .catch((error) => console.log(error));
  }, [cellTypeId]);

  const tableRows: TableRow[] = [];
  for (const collection of collections) {
    const tissueNames = collection.tissue.map((tissue) => tissue.label);
    const diseaseNames = collection.disease.map((disease) => disease.label);
    const organismNames = collection.organism.map((organism) => organism.label);

    tableRows.push({
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

  return (
    <>
      <TableTitleWrapper id={SOURCE_DATA_SECTION_ID}>
        <TableTitle>Data</TableTitle>
      </TableTitleWrapper>
      {tableRows.length && (
        <Table<TableRow> columns={tableColumns} rows={tableRows} />
      )}
    </>
  );
};
export default SourceDataTable;
