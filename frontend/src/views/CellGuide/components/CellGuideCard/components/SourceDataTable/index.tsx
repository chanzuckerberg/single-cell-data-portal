import React, { ReactElement, useMemo } from "react";
import { Tooltip } from "@czi-sds/components";
import {
  TableTitle,
  TableTitleWrapper,
  TableUnavailableContainer,
  TableUnavailableHeader,
  TableUnavailableDescription,
} from "../common/style";
import Table from "../common/Table";
import Link from "../common/Link";
import { StyledTag } from "./style";
import { useSourceData } from "src/common/queries/cellGuide";

export const CELL_GUIDE_CARD_SOURCE_DATA_TABLE =
  "cell-guide-card-source-data-table";

interface TableRow {
  collection: ReactElement;
  publication: ReactElement | string;
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
    const sortedCollections = collections.sort((a, b) => {
      const aOrganisms = a.organism.map((organism) => organism.label);
      const bOrganisms = b.organism.map((organism) => organism.label);
      if (aOrganisms.length === 0 && bOrganisms.length === 0) return 0;
      if (aOrganisms.includes("Homo sapiens")) return -1;
      if (bOrganisms.includes("Homo sapiens")) return 1;
      const aFirstOrganism = aOrganisms.at(0) ?? "";
      const bFirstOrganism = bOrganisms.at(0) ?? "";
      return aFirstOrganism.localeCompare(bFirstOrganism);
    });
    for (const collection of sortedCollections) {
      const tissueNames = collection.tissue.map((tissue) => tissue.label);
      const diseaseNames = collection.disease.map((disease) => disease.label);
      const organismNames = collection.organism.map(
        (organism) => organism.label
      );
      rows.push({
        collection: (
          <Link
            key={`collection-name-${collection.collection_name}-${index}`}
            label={collection.collection_name}
            url={collection.collection_url}
          />
        ),
        publication: collection.publication_url ? (
          <Link
            key={`publication-url-${collection.publication_title}-${index}`}
            label={collection.publication_title}
            url={`https://doi.org/${collection.publication_url}`}
          />
        ) : (
          "No publication"
        ),
        tissue:
          tissueNames.length <= 2 ? (
            <div>
              {tissueNames.map((tissue) => {
                return <div key={`tissue-${tissue}-${index}`}>{tissue}</div>;
              })}
            </div>
          ) : (
            <Tooltip
              sdsStyle="dark"
              placement="right"
              width="default"
              arrow
              title={
                <div>
                  {tissueNames.map((tissue) => {
                    return (
                      <div key={`tissue-${tissue}-${index}`}>{tissue}</div>
                    );
                  })}
                </div>
              }
            >
              <span>
                <StyledTag
                  color="gray"
                  sdsType="secondary"
                  label={`${tissueNames.length} tissues`}
                />
              </span>
            </Tooltip>
          ),
        disease: (
          <div>
            {diseaseNames.map((disease) => {
              return <div key={`disease-${disease}-${index}`}>{disease}</div>;
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
    <div data-testid={CELL_GUIDE_CARD_SOURCE_DATA_TABLE}>
      <TableTitleWrapper>
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
    </div>
  );
};
export default SourceDataTable;
