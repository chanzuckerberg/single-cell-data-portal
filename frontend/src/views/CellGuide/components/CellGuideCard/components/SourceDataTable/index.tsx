import React, { ReactElement, useEffect, useMemo, useState } from "react";
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
import { Pagination } from "@mui/material";

import { CELL_GUIDE_CARD_SOURCE_DATA_TABLE } from "src/views/CellGuide/components/CellGuideCard/components/SourceDataTable/constants";

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

const ROWS_PER_PAGE = 10;

const SourceDataTable = ({ cellTypeId }: Props) => {
  const { data: collections } = useSourceData(cellTypeId);
  const [page, setPage] = useState(1);

  const handlePageChange = (
    _event: React.ChangeEvent<unknown>,
    page: number
  ) => {
    setPage(page);
  };

  // Handle cell type change, set marker genes table page back to 1
  useEffect(() => {
    setPage(1);
  }, [cellTypeId]);

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
              sdsStyle="light"
              placement="top"
              width="wide"
              leaveDelay={0}
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
          <>
            <div>
              {diseaseNames.length <= 2 ? (
                <div>
                  {diseaseNames.map((disease) => {
                    return (
                      <div key={`disease-${disease}-${index}`}>{disease}</div>
                    );
                  })}
                </div>
              ) : (
                <>
                  {/* If 'normal' exists then have it outside of the overflow tag */}
                  {diseaseNames.includes("normal") && (
                    <div key={`disease-normal`}>normal</div>
                  )}

                  <Tooltip
                    sdsStyle="light"
                    placement="top"
                    width="wide"
                    leaveDelay={0}
                    title={
                      <div>
                        {diseaseNames
                          .filter((disease) => disease !== "normal")
                          .map((disease) => {
                            return (
                              <div key={`disease-${disease}-${index}`}>
                                {disease}
                              </div>
                            );
                          })}
                      </div>
                    }
                  >
                    <span>
                      <StyledTag
                        color="gray"
                        sdsType="secondary"
                        label={`${
                          diseaseNames.includes("normal")
                            ? diseaseNames.length - 1
                            : diseaseNames.length
                        } diseases`}
                      />
                    </span>
                  </Tooltip>
                </>
              )}
            </div>
          </>
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

  const pageCount = Math.ceil(tableRows.length / ROWS_PER_PAGE);

  const tableUnavailableComponent = (
    <TableUnavailableContainer>
      <TableUnavailableHeader>Source data unavailable</TableUnavailableHeader>
      <TableUnavailableDescription>
        Neither this cell type nor its ancestors and descendants are available
        in any collections.
      </TableUnavailableDescription>
    </TableUnavailableContainer>
  );

  const tableComponent = (
    <Table<TableRow>
      columns={tableColumns}
      rows={
        tableRows.slice(
          (page - 1) * ROWS_PER_PAGE,
          page * ROWS_PER_PAGE
        ) as TableRow[]
      }
    />
  );

  return (
    <div data-testid={CELL_GUIDE_CARD_SOURCE_DATA_TABLE}>
      <TableTitleWrapper>
        <TableTitle>Data</TableTitle>
      </TableTitleWrapper>

      {tableRows.length > 0 ? (
        <div>
          {tableComponent}
          <Pagination
            count={pageCount}
            page={page}
            onChange={handlePageChange}
          />
        </div>
      ) : (
        tableUnavailableComponent
      )}
    </div>
  );
};
export default SourceDataTable;
