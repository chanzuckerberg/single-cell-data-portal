import React, {
  ReactElement,
  useEffect,
  useMemo,
  useState,
  SetStateAction,
  Dispatch,
} from "react";
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
import { useIsComponentPastBreakpointWidth } from "../common/hooks/useIsComponentPastBreakpoint";
import {
  CELL_GUIDE_CARD_SOURCE_DATA_TABLE,
  SOURCE_DATA_TABLE_BREAKPOINT_PX,
} from "./constants";
import { useDataSourceFilter } from "./hooks/useDataSourceFilter";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";

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
  organId: string;
  organismName: string;
  skinnyMode: boolean;
  setTooltipContent: Dispatch<
    SetStateAction<{
      title: string;
      element: JSX.Element;
    } | null>
  >;
}

const ROWS_PER_PAGE = 10;

const SourceDataTable = ({
  cellTypeId,
  organId,
  organismName,
  skinnyMode,
  setTooltipContent,
}: Props) => {
  const { data: collections } = useSourceData(cellTypeId);
  const [page, setPage] = useState(1);
  const { isPastBreakpoint, containerRef } = useIsComponentPastBreakpointWidth(
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
    selectedOrganId: organId,
  });

  const tableRows: TableRow[] = useMemo(() => {
    const createTooltipClickHandler = (names: string[], title: string) => {
      return () => {
        if (skinnyMode) {
          setTooltipContent({
            title,
            element: (
              <div>
                {names.map((name, index) => (
                  <div key={`name-${name}-${index}`}>{name}</div>
                ))}
              </div>
            ),
          });
        }
      };
    };
    return filteredCollections.map((collections, index) =>
      createTableRow(
        collections,
        index,
        isPastBreakpoint,
        skinnyMode,
        createTooltipClickHandler
      )
    );
  }, [filteredCollections, isPastBreakpoint, skinnyMode, setTooltipContent]);

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
      {tableRows
        .slice((page - 1) * ROWS_PER_PAGE, page * ROWS_PER_PAGE)
        .map((row, index) => {
          return (
            <MobileSourceDataTableEntry key={index} highlight={index % 2 === 1}>
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

function createTissueContent(
  tissueNames: string[],
  isPastBreakpoint: boolean,
  skinnyMode: boolean,
  createTooltipClickHandler: (names: string[], title: string) => () => void
) {
  return tissueNames.length <= (isPastBreakpoint ? 1 : 3) ? (
    tissueNames.map((tissue, index) => (
      <span key={`tissue-${tissue}-${index}`}>
        <StyledTag color="gray" sdsType="secondary" label={tissue} />
      </span>
    ))
  ) : (
    <Tooltip
      sdsStyle="light"
      placement="top"
      width="wide"
      disableHoverListener={skinnyMode}
      leaveDelay={0}
      title={tissueNames.map((tissue, index) => (
        <div key={`tissue-${tissue}-${index}`}>{tissue}</div>
      ))}
    >
      <span onClick={createTooltipClickHandler(tissueNames, "Tissues")}>
        <StyledTag
          color="gray"
          sdsType="secondary"
          label={generateTagLabel(tissueNames, "tissue", "tissues")}
        />
      </span>
    </Tooltip>
  );
}

function createDiseaseContent(
  diseaseNames: string[],
  isPastBreakpoint: boolean,
  skinnyMode: boolean,
  createTooltipClickHandler: (names: string[], title: string) => () => void
) {
  return diseaseNames.length <= (isPastBreakpoint ? 1 : 3) ? (
    diseaseNames
      .sort((a, b) => (a === "normal" ? -1 : b === "normal" ? 1 : 0))
      .map((disease, index) => (
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
        disableHoverListener={skinnyMode}
        title={diseaseNames
          .filter((disease) => isPastBreakpoint || disease !== "normal")
          .map((disease, index) => (
            <div key={`disease-${disease}-${index}`}>{disease}</div>
          ))}
      >
        <span onClick={createTooltipClickHandler(diseaseNames, "Diseases")}>
          <StyledTag
            color="gray"
            sdsType="secondary"
            label={generateTagLabel(diseaseNames, "disease", "diseases")}
          />
        </span>
      </Tooltip>
    </>
  );
}

const createOrganismContent = (
  organismNames: string[],
  isPastBreakpoint: boolean,
  skinnyMode: boolean,
  createTooltipClickHandler: (names: string[], title: string) => () => void
) => {
  return organismNames.length <= (isPastBreakpoint ? 1 : 3) ? (
    organismNames.map((organism, index) => (
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
        disableHoverListener={skinnyMode}
        title={organismNames.map((organism, index) => (
          <div key={`organism-${organism}-${index}`}>{organism}</div>
        ))}
      >
        <span onClick={createTooltipClickHandler(organismNames, "Organisms")}>
          <StyledTag
            color="gray"
            sdsType="secondary"
            label={generateTagLabel(organismNames, "organism", "organisms")}
          />
        </span>
      </Tooltip>
    </>
  );
};

function createTableRow(
  collection: SourceCollectionsQueryResponseEntry,
  index: number,
  isPastBreakpoint: boolean,
  skinnyMode: boolean,
  createTooltipClickHandler: (names: string[], title: string) => () => void
) {
  const tissueNames = collection.tissue.map((tissue) => tissue.label);
  const diseaseNames = collection.disease.map((disease) => disease.label);
  const organismNames = collection.organism.map((organism) => organism.label);

  const tissueContent = createTissueContent(
    tissueNames,
    isPastBreakpoint,
    skinnyMode,
    createTooltipClickHandler
  );

  const diseaseContent = createDiseaseContent(
    diseaseNames,
    isPastBreakpoint,
    skinnyMode,
    createTooltipClickHandler
  );

  const organismContent = createOrganismContent(
    organismNames,
    isPastBreakpoint,
    skinnyMode,
    createTooltipClickHandler
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
  return (
    <Link
      key={key}
      label={label}
      url={url}
      onClick={() =>
        track(EVENTS.CG_DATA_COLLECTION_CLICKED, { collection: label })
      }
    />
  );
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
