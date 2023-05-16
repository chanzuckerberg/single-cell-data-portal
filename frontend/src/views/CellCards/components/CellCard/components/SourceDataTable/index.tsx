import React, { ReactElement, useMemo } from "react";
import {
  TableTitle,
  TableTitleWrapper,
  TableUnavailableContainer,
  TableUnavailableHeader,
  TableUnavailableDescription,
} from "../common/style";
import Table from "../Table";
import Link from "../common/Link";
import { StyledTag } from "./style";
import { useSourceData } from "src/common/queries/cellCards";
import { SOURCE_DATA_SECTION_ID } from "../CellCardSidebar";

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
    let index = 0;
    for (const collection of collections) {
      const tissueNames = collection.tissue.map((tissue) => tissue.label);
      const diseaseNames = collection.disease.map((disease) => disease.label);
      const organismNames = collection.organism.map(
        (organism) => organism.label
      );

      rows.push({
        collection: (
          <Link
            key={`collection-name-${collection.collection_name}-${index}`}
            title={collection.collection_name}
            url={collection.collection_url}
          />
        ),
        publication: (
          <Link
            key={`publication-url-${collection.publication_title}-${index}`}
            title={collection.publication_title}
            url={`https://doi.org/${collection.publication_url}`}
          />
        ),
        tissue:
          tissueNames.length <= 7 ? (
            <div>
              {tissueNames.map((tissue) => {
                return (
                  <div key={`tissue-${tissue}-${index}`}>
                    {tissue.charAt(0).toUpperCase() + tissue.slice(1)}
                  </div>
                );
              })}
            </div>
          ) : (
            <StyledTag
              color="gray"
              sdsType="secondary"
              label={`${tissueNames.length} tissues`}
            />
          ),
        disease: (
          <div>
            {diseaseNames.map((disease) => {
              return (
                <div key={`disease-${disease}-${index}`}>
                  {disease.charAt(0).toUpperCase() + disease.slice(1)}
                </div>
              );
            })}
          </div>
        ),
        organism: (
          <div>
            {organismNames.map((organism) => {
              return (
                <div key={`organism-${organism}-${index}`}>{organism}</div>
              );
            })}
          </div>
        ),
      });
      index = index + 1;
    }
    return rows;
  }, [collections]);

  return (
    <>
      <TableTitleWrapper id={SOURCE_DATA_SECTION_ID}>
        <TableTitle>Data</TableTitle>
      </TableTitleWrapper>
      {tableRows.length > 0 ? (
        <Table<TableRow> columns={tableColumns} rows={tableRows} />
      ) : (
        <TableUnavailableContainer>
          <TableUnavailableHeader>
            Source data unavailable
          </TableUnavailableHeader>
          <TableUnavailableDescription>
            Neither this cell type nor its ancestors and descendants are
            available in any collections.
          </TableUnavailableDescription>
        </TableUnavailableContainer>
      )}
    </>
  );
};
export default SourceDataTable;
