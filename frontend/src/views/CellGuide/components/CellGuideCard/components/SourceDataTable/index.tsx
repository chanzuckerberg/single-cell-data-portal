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
import {
  SourceDataTableWrapper,
  StyledTag,
  MobileSourceDataTableWrapper,
  MobileSourceDataTableEntry,
  MobileSourceDataTableEntryRow,
} from "./style";
import {
  SourceCollectionsQueryResponseEntry,
  useSourceData,
} from "src/common/queries/cellGuide";
import { Pagination } from "@mui/material";
import useIsComponentPastBreakpoint from "../common/hooks/useIsComponentPastBreakpoint";
import {
  CELL_GUIDE_CARD_SOURCE_DATA_TABLE,
  SOURCE_DATA_TABLE_BREAKPOINT_PX,
} from "./constants";
import { useDataSourceFilter } from "./hooks/useDataSourceFilter";

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
  organName: string;
  organId: string;
  organismName: string;
}

const ROWS_PER_PAGE = 10;

const SourceDataTable = ({
  cellTypeId,
  organName,
  organId,
  organismName,
}: Props) => {
  const { data: collections } = useSourceData(cellTypeId);
  const [page, setPage] = useState(1);
  const { isPastBreakpoint, containerRef } = useIsComponentPastBreakpoint(
    SOURCE_DATA_TABLE_BREAKPOINT_PX
  );

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

  const filteredCollections = useDataSourceFilter({
    collections: collections ?? [],
    selectedOrganismLabel: organismName,
    selectedOrganLabel: organName,
    selectedOrganId: organId,
  });

  const tableRows: TableRow[] = useMemo(() => {
    return filteredCollections.map((collections, index) =>
      createTableRow(collections, index, isPastBreakpoint)
    );
  }, [filteredCollections, isPastBreakpoint]);

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

  const tableComponent = isPastBreakpoint ? (
    <MobileSourceDataTableWrapper>
      {tableRows.map((row, index) => {
        return (
          <MobileSourceDataTableEntry key={index} index={index}>
            {row.collection}
            {row.publication}
            <MobileSourceDataTableEntryRow>
              {row.tissue}
              {row.disease}
              {row.organism}
            </MobileSourceDataTableEntryRow>
          </MobileSourceDataTableEntry>
        );
      })}
    </MobileSourceDataTableWrapper>
  ) : (
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
    <SourceDataTableWrapper data-testid={CELL_GUIDE_CARD_SOURCE_DATA_TABLE}>
      <TableTitleWrapper>
        <TableTitle>Data</TableTitle>
      </TableTitleWrapper>

      {tableRows.length > 0 ? (
        <div ref={containerRef}>
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
    </SourceDataTableWrapper>
  );
};

function createTableRow(
  collection: SourceCollectionsQueryResponseEntry,
  index: number,
  isPastBreakpoint: boolean
) {
  const tissueNames = collection.tissue.map((tissue) => tissue.label);
  const diseaseNames = collection.disease.map((disease) => disease.label);
  const organismNames = collection.organism.map((organism) => organism.label);

  const tissueContent =
    tissueNames.length <= (isPastBreakpoint ? 1 : 3) ? (
      tissueNames.map((tissue) => (
        <span key={`tissue-${tissue}-${index}`}>
          <StyledTag color="gray" sdsType="secondary" label={tissue} />
        </span>
      ))
    ) : (
      <Tooltip
        sdsStyle="light"
        placement="top"
        width="wide"
        leaveDelay={0}
        title={tissueNames.map((tissue) => (
          <div key={`tissue-${tissue}-${index}`}>{tissue}</div>
        ))}
      >
        <span>
          <StyledTag
            color="gray"
            sdsType="secondary"
            label={generateTagLabel(tissueNames, "tissue", "tissues")}
          />
        </span>
      </Tooltip>
    );

  const diseaseContent =
    diseaseNames.length <= (isPastBreakpoint ? 1 : 3) ? (
      diseaseNames
        .sort((a, b) => (a === "normal" ? -1 : b === "normal" ? 1 : 0))
        .map((disease) => (
          <span key={`disease-${disease}-${index}`}>
            {<StyledTag color="gray" sdsType="secondary" label={disease} />}
          </span>
        ))
    ) : (
      <>
        {/* If 'normal' exists then have it outside of the overflow tag */}
        {diseaseNames.includes("normal") && !isPastBreakpoint && (
          <span>
            <StyledTag color="gray" sdsType="secondary" label="normal" />
          </span>
        )}
        <Tooltip
          sdsStyle="light"
          placement="top"
          width="wide"
          leaveDelay={0}
          title={diseaseNames
            .filter((disease) => isPastBreakpoint || disease !== "normal")
            .map((disease) => (
              <div key={`disease-${disease}-${index}`}>{disease}</div>
            ))}
        >
          <span>
            <StyledTag
              color="gray"
              sdsType="secondary"
              label={generateTagLabel(diseaseNames, "disease", "diseases")}
            />
          </span>
        </Tooltip>
      </>
    );

  const organismContent =
    organismNames.length <= (isPastBreakpoint ? 1 : 3) ? (
      organismNames.map((organism) => (
        <span key={`organism-${organism}-${index}`}>
          {<StyledTag color="gray" sdsType="secondary" label={organism} />}
        </span>
      ))
    ) : (
      <>
        <Tooltip
          sdsStyle="light"
          placement="top"
          width="wide"
          leaveDelay={0}
          title={organismNames.map((organism) => (
            <div key={`organism-${organism}-${index}`}>{organism}</div>
          ))}
        >
          <span>
            <StyledTag
              color="gray"
              sdsType="secondary"
              label={generateTagLabel(organismNames, "organism", "organisms")}
            />
          </span>
        </Tooltip>
      </>
    );

  return {
    collection: generateLink(
      `collection-name-${collection.collection_name}-${index}`,
      collection.collection_name,
      collection.collection_url
    ),
    publication: <div>{collection.publication_title}</div>,
    tissue: <div>{tissueContent}</div>,
    disease: <div>{diseaseContent}</div>,
    organism: <div>{organismContent}</div>,
  };
}

function generateLink(key: string, label: string, url: string) {
  return <Link key={key} label={label} url={url} />;
}

function generateTagLabel(
  names: string[],
  singularLabel: string,
  pluralLabel: string
) {
  if (names.includes(singularLabel)) {
    return singularLabel;
  } else {
    return `${names.length} ${pluralLabel}`;
  }
}

export default SourceDataTable;
