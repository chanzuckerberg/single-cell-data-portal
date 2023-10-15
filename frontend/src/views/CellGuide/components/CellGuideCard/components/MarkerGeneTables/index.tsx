/* eslint-disable react/no-unescaped-entities */
import React, {
  Dispatch,
  ReactElement,
  ReactNode,
  SetStateAction,
  useCallback,
  useEffect,
  useMemo,
  useState,
} from "react";
import { ButtonIcon, Tooltip } from "@czi-sds/components";
import {
  TableTitle,
  TableTitleWrapper,
  TableUnavailableContainer,
  TableUnavailableHeader,
  TableUnavailableDescription,
  TableTitleInnerWrapper,
  FlexRow,
} from "../common/style";
import Link from "../common/Link";
import {
  PublicationLinkWrapper,
  TableSelectorButton,
  TableSelectorRow,
  TableTitleOuterWrapper,
  StyledHeadCellContent,
  MarkerStrengthContainer,
  MarkerGenePagination,
  MarkerGeneInfo,
  MarkerGeneTooltipText,
  MarkerGeneTooltipSubtext,
  MarkerGeneTableWrapper,
  StyledCellNumerical,
  NoWrapWrapper,
  StyledLink,
  ReferenceTooltipWrapper,
  StyledImageWrapper,
} from "./style";
import Table from "../common/Table";
import { Pagination } from "@mui/material";
import {
  CanonicalMarkersQueryResponse,
  ComputationalMarkersQueryResponse,
  useAllTissuesLookupTables,
  useCanonicalMarkers,
  useComputationalMarkers,
} from "src/common/queries/cellGuide";
import {
  ComputationalMarkerGeneTableData,
  useComputationalMarkerGenesTableRowsAndFilters,
} from "./hooks/computational_markers";
import treeDendrogram from "src/common/images/TreeDendogram.svg";
import {
  CanonicalMarkerGeneTableData,
  useCanonicalMarkerGenesTableRowsAndFilters,
} from "./hooks/canonical_markers";
import { useIsComponentPastBreakpointWidth } from "../common/hooks/useIsComponentPastBreakpoint";
import HelpTooltip from "../common/HelpTooltip";
import { ROUTES } from "src/common/constants/routes";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { CENSUS_LINK } from "src/components/Header/components/Nav";
import {
  CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE,
  CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE_SELECTOR,
  CELL_GUIDE_CARD_ENRICHED_GENES_TABLE,
  CELL_GUIDE_CARD_ENRICHED_GENES_TABLE_SELECTOR,
  EXPRESSION_SCORE_TOOLTIP_TEST_ID,
  MARKER_GENES_CANONICAL_TOOLTIP_TEST_ID,
  MARKER_GENES_COMPUTATIONAL_TOOLTIP_TEST_ID,
  MARKER_SCORE_TOOLTIP_TEST_ID,
  PERCENT_OF_CELLS_TOOLTIP_TEST_ID,
  MARKER_GENES_CANONICAL_BREAKPOINT_PX,
  MARKER_GENES_COMPUTATIONAL_BREAKPOINT_PX,
  MARKER_GENES_TREE_ICON_BUTTON_TEST_ID,
} from "src/views/CellGuide/components/CellGuideCard/components/MarkerGeneTables/constants";
import { FMG_GENE_STRENGTH_THRESHOLD } from "src/views/WheresMyGene/common/constants";
import Image from "next/image";
import { CellType } from "../../../common/OntologyDagView/common/types";

function getEmptyComputationalMarkerGenesTableUIMessageDetail(
  allFilteredByLowMarkerScore: boolean
): string {
  if (allFilteredByLowMarkerScore) {
    return (
      "There are no computational marker genes with marker scores above the threshold of " +
      FMG_GENE_STRENGTH_THRESHOLD +
      "."
    );
  }

  return "Computational marker genes for this cell type are unavailable at this time.";
}

function getEmptyCanonicalMarkerGenesTableUIMessageDetail(): string {
  return "Canonical marker genes for this cell type are unavailable at this time.";
}

const ROWS_PER_PAGE = 10;

// Computational marker gene table types
interface TableRowEnrichedGenes {
  symbol: ReactNode;
  name: ReactNode;
  marker_score: ReactNode;
  me: ReactNode;
  pc: ReactNode;
}

// Canonical marker gene table types
interface TableRowCanonicalGenes {
  symbol: ReactNode;
  name: string;
  references: ReactNode;
}

// Table row type
type TableRow = (TableRowEnrichedGenes | TableRowCanonicalGenes) & {
  symbolId: string;
};

// Tooltip components
const enrichedGenesTooltipComponent = (
  <div>
    <MarkerGeneTooltipText>
      {"Computational marker genes are derived from the "}
      <Link label={"CELLxGENE Census"} url={CENSUS_LINK} />
      {". They are computed utilizing the same methodology as featured in our "}
      <Link
        label={"Find Marker Genes feature from the Gene Expression application"}
        url={ROUTES.FMG_DOCS}
      />
      {"."}
    </MarkerGeneTooltipText>
  </div>
);
const canonicalMarkerGenesTooltipComponent = (
  <div>
    <MarkerGeneTooltipText>
      {
        "Canonical marker genes and associated publications were derived from the "
      }
      <Link
        label={"Anatomical Structures, Cell Types and Biomarkers (ASCT+B)"}
        url={"https://humanatlas.io/asctb-tables"}
      />
      {
        " tables from the 5th Human Reference Atlas release (July 2023). The tables are authored and reviewed by an international team of anatomists, pathologists, physicians, and other experts."
      }
    </MarkerGeneTooltipText>
    <br />
    <MarkerGeneTooltipSubtext>
      <Link
        label={
          <i>
            Börner, Katy, et al. "Anatomical structures, cell types and
            biomarkers of the Human Reference Atlas." Nature cell biology 23.11
            (2021): 1117-1128.
          </i>
        }
        url={"https://www.nature.com/articles/s41556-021-00788-6"}
      />
    </MarkerGeneTooltipSubtext>
  </div>
);

interface Props {
  cellTypeId: string;
  setGeneInfoGene: React.Dispatch<React.SetStateAction<string | null>>;
  cellTypeName: string;
  skinnyMode: boolean;
  setTooltipContent: Dispatch<
    SetStateAction<{
      title: string;
      element: JSX.Element;
    } | null>
  >;
  organName: string;
  organId: string;
  organismName: string;
  selectedGene?: string;
  selectGene: (gene: string) => void;
  setCellInfoCellType?: React.Dispatch<React.SetStateAction<CellType | null>>;
  cellInfoCellType?: CellType | null;
}

const MarkerGeneTables = ({
  cellTypeId,
  cellTypeName,
  setGeneInfoGene,
  skinnyMode,
  setTooltipContent,
  organName,
  organId,
  organismName,
  selectedGene,
  selectGene,
  cellInfoCellType,
  setCellInfoCellType,
}: Props) => {
  // 0 is canonical marker genes, 1 is computational marker genes
  const [activeTable, setActiveTable] = useState(1);
  const [computationalMarkerGenes, setComputationalMarkerGenes] =
    useState<ComputationalMarkersQueryResponse>([]);

  const { isPastBreakpoint, containerRef } = useIsComponentPastBreakpointWidth(
    activeTable
      ? MARKER_GENES_COMPUTATIONAL_BREAKPOINT_PX
      : MARKER_GENES_CANONICAL_BREAKPOINT_PX
  );

  const [canonicalMarkerGenes, setCanonicalMarkerGenes] =
    useState<CanonicalMarkersQueryResponse>([]);

  const { data: enrichedGenes } = useComputationalMarkers(cellTypeId);
  const { data: canonicalMarkers } = useCanonicalMarkers(cellTypeId);

  // Computational marker gene table column names
  const tableColumnNamesEnrichedGenes: Record<
    keyof TableRowEnrichedGenes,
    ReactElement | string
  > = useMemo(
    () => ({
      symbol: "Symbol",
      name: "Name",
      marker_score: (
        <div>
          <StyledHeadCellContent>
            Marker Score
            <HelpTooltip
              skinnyMode={skinnyMode}
              title="Marker Score"
              setTooltipContent={setTooltipContent}
              dark
              buttonDataTestId={MARKER_SCORE_TOOLTIP_TEST_ID}
              text={
                <>
                  Marker score interpretation:
                  <br />
                  <MarkerStrengthContainer>
                    {"Low: <1 | Medium: 1-2 | High: >2"}
                  </MarkerStrengthContainer>
                  <br />
                  <div>
                    Marker genes are highly and uniquely expressed in the cell
                    type relative to all other cell types.
                  </div>
                  <br />
                  <div>
                    <a href={ROUTES.FMG_DOCS} rel="noopener" target="_blank">
                      Click to read more about the identification method.
                    </a>
                  </div>
                </>
              }
            />
          </StyledHeadCellContent>
        </div>
      ),
      me: (
        <StyledHeadCellContent>
          {!isPastBreakpoint ? "Expression Score" : "Exp. Score"}
          <HelpTooltip
            skinnyMode={skinnyMode}
            title="Expression Score"
            setTooltipContent={setTooltipContent}
            dark
            buttonDataTestId={EXPRESSION_SCORE_TOOLTIP_TEST_ID}
            text={
              <div>
                The expression score is the average{" "}
                <a
                  href={ROUTES.WMG_DOCS_DATA_PROCESSING}
                  target="_blank"
                  rel="noreferrer noopener"
                >
                  rankit-normalized gene expression
                </a>{" "}
                among cells in the cell type that have non-zero values.
              </div>
            }
          />
        </StyledHeadCellContent>
      ),
      pc: (
        <StyledHeadCellContent>
          % of Cells
          <HelpTooltip
            skinnyMode={skinnyMode}
            title="% of Cells"
            setTooltipContent={setTooltipContent}
            dark
            buttonDataTestId={PERCENT_OF_CELLS_TOOLTIP_TEST_ID}
            text={
              <div>
                Percentage of cells expressing a gene in the cell type. These
                numbers are calculated after cells with{" "}
                <a
                  href={ROUTES.WMG_DOCS_DATA_PROCESSING}
                  target="_blank"
                  rel="noreferrer noopener"
                >
                  low coverage and low expression values
                </a>{" "}
                have been filtered out.
              </div>
            }
          />
        </StyledHeadCellContent>
      ),
    }),
    [setTooltipContent, skinnyMode, isPastBreakpoint]
  );

  const allTissuesLabelToIdMap = useAllTissuesLookupTables(cellTypeId);

  useEffect(() => {
    if (enrichedGenes) {
      setComputationalMarkerGenes(enrichedGenes);
    }
  }, [enrichedGenes]);

  useEffect(() => {
    if (canonicalMarkers) {
      setCanonicalMarkerGenes(canonicalMarkers);
    }
  }, [canonicalMarkers]);

  const [page, setPage] = useState(1);

  const { computationalMarkerGeneTableData, allFilteredByLowMarkerScore } =
    useComputationalMarkerGenesTableRowsAndFilters({
      genes: computationalMarkerGenes,
      allTissuesLabelToIdMap: allTissuesLabelToIdMap,
      selectedOrganId: organId,
      selectedOrganismLabel: organismName,
    });

  const { canonicalMarkerGeneTableData } =
    useCanonicalMarkerGenesTableRowsAndFilters({
      genes: canonicalMarkerGenes,
      allTissuesLabelToIdMap: allTissuesLabelToIdMap,
      selectedOrganLabel: organName,
      selectedOrganId: organId,
      selectedOrganismLabel: organismName,
    });

  const getSymbol = useCallback(
    (
      row: ComputationalMarkerGeneTableData | CanonicalMarkerGeneTableData,
      showEye = false
    ) => (
      <NoWrapWrapper isSelected={row.symbol === selectedGene}>
        {row.symbol}{" "}
        <ButtonIcon
          aria-label={`display gene info for ${row.symbol}`}
          className="hover-button"
          sdsIcon="infoCircle"
          sdsSize="small"
          sdsType="secondary"
          onClick={() => setGeneInfoGene(row.symbol.toUpperCase())}
        />
        {showEye && (
          <StyledImageWrapper
            className="hover-button"
            isActive={row.symbol === selectedGene}
            onClick={() => {
              skinnyMode && setCellInfoCellType && setCellInfoCellType(null);
              track(EVENTS.CG_MARKER_GENE_MODE_CLICKED, {
                cell_type: setCellInfoCellType
                  ? cellInfoCellType?.cellTypeName
                  : cellTypeName,
                inSideBar: !!cellInfoCellType,
              });
              selectGene(row.symbol);
            }}
          >
            <Image
              data-testid={MARKER_GENES_TREE_ICON_BUTTON_TEST_ID(row.symbol)}
              src={treeDendrogram}
              alt={`activate marker gene mode for ${row.symbol}`}
              width="12px"
              height="12px"
            />
          </StyledImageWrapper>
        )}
      </NoWrapWrapper>
    ),
    [
      selectedGene,
      setGeneInfoGene,
      selectGene,
      setCellInfoCellType,
      skinnyMode,
      cellInfoCellType,
      cellTypeName,
    ]
  );

  const tableRows: TableRow[] = useMemo(() => {
    const referenceClickHandlerMobileView = (
      row: CanonicalMarkerGeneTableData
    ) => {
      return (event: React.MouseEvent<HTMLSpanElement>) => {
        event.preventDefault();
        setTooltipContent({
          title: `${row.symbol} marker gene references`,
          element: (
            <ReferenceTooltipWrapper>
              {row.referenceData.publicationTitles.map(
                (publicationTitle, index) => {
                  const keyVal = `${row.referenceData.publications[index]}-${index}}`;
                  const referenceIndexLabel =
                    (row.referenceData.publicationTitlesToIndex.get(
                      publicationTitle
                    ) ?? 0) + 1;
                  return (
                    <div key={`${keyVal}-tooltip`}>
                      <StyledLink
                        key={`${keyVal}`}
                        label={`[${referenceIndexLabel}] ${publicationTitle
                          .split("\n\n")
                          .at(0)}`}
                        url={`https://doi.org/${row.referenceData.publications[index]}`}
                      />
                      <br />
                      <i>{publicationTitle.split("\n\n").at(1)}</i>
                    </div>
                  );
                }
              )}
            </ReferenceTooltipWrapper>
          ),
        });
      };
    };

    return activeTable
      ? computationalMarkerGeneTableData.map((row) => ({
          ...row,
          me: <StyledCellNumerical> {row.me} </StyledCellNumerical>,
          pc: <StyledCellNumerical> {row.pc} </StyledCellNumerical>,
          marker_score: (
            <StyledCellNumerical> {row.marker_score} </StyledCellNumerical>
          ),
          symbolId: row.symbol,
          symbol: getSymbol(row, true),
        }))
      : canonicalMarkerGeneTableData.map((row) => ({
          ...row,
          symbolId: row.symbol,
          symbol: getSymbol(row),
          references: (
            <PublicationLinkWrapper>
              {row.referenceData.publicationTitles.map(
                (publicationTitle, index) => {
                  if (
                    publicationTitle &&
                    row.referenceData.publications[index]
                  ) {
                    const referenceIndexLabel =
                      (row.referenceData.publicationTitlesToIndex.get(
                        publicationTitle
                      ) ?? 0) + 1;
                    const keyVal = `${row.referenceData.publications[index]}-${index}}`;
                    return (
                      <Tooltip
                        key={`${keyVal}-tooltip`}
                        placement="top"
                        disableHoverListener={skinnyMode}
                        width="default"
                        arrow={false}
                        title={
                          <div>
                            {publicationTitle.split("\n\n").at(0)}
                            <br />
                            <br />
                            <i>{publicationTitle.split("\n\n").at(1)}</i>
                          </div>
                        }
                        leaveDelay={0}
                      >
                        <span
                          key={`${keyVal}-span`}
                          onClick={
                            skinnyMode
                              ? referenceClickHandlerMobileView(row)
                              : undefined
                          }
                        >
                          <StyledLink
                            key={`${keyVal}`}
                            label={`[${referenceIndexLabel}]`}
                            url={`https://doi.org/${row.referenceData.publications[index]}`}
                          />
                        </span>
                      </Tooltip>
                    );
                  }
                }
              )}
            </PublicationLinkWrapper>
          ),
        }));
  }, [
    activeTable,
    canonicalMarkerGeneTableData,
    computationalMarkerGeneTableData,
    setTooltipContent,
    getSymbol,
    skinnyMode,
  ]);

  const genesForShareUrl = `${tableRows
    .map((row) => row.symbolId)
    .join("%2C")}&cellTypes=${cellTypeName.replace(" ", "+")}`;

  const pageCount = Math.ceil(tableRows.length / ROWS_PER_PAGE);
  const tableComponent = useMemo(() => {
    const tableColumnsEnrichedGenes: Array<keyof TableRowEnrichedGenes> =
      !isPastBreakpoint
        ? ["symbol", "name", "marker_score", "me", "pc"]
        : ["symbol", "marker_score", "me", "pc"];

    const tableColumnsCanonicalGenes: Array<keyof TableRowCanonicalGenes> =
      !isPastBreakpoint
        ? ["symbol", "name", "references"]
        : ["symbol", "references"];
    // Canonical marker gene table column names
    const tableColumnNamesCanonicalGenes: Record<
      keyof TableRowCanonicalGenes,
      string
    > = {
      symbol: "Symbol",
      name: "Name",
      references: "References",
    };
    return activeTable ? (
      <Table<TableRowEnrichedGenes>
        testId={CELL_GUIDE_CARD_ENRICHED_GENES_TABLE}
        columns={tableColumnsEnrichedGenes}
        rows={
          tableRows.slice(
            (page - 1) * ROWS_PER_PAGE,
            page * ROWS_PER_PAGE
          ) as TableRowEnrichedGenes[]
        }
        columnIdToName={tableColumnNamesEnrichedGenes}
      />
    ) : (
      <Table<TableRowCanonicalGenes>
        testId={CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE}
        columns={tableColumnsCanonicalGenes}
        rows={
          tableRows.slice(
            (page - 1) * ROWS_PER_PAGE,
            page * ROWS_PER_PAGE
          ) as TableRowCanonicalGenes[]
        }
        columnIdToName={tableColumnNamesCanonicalGenes}
      />
    );
  }, [
    activeTable,
    page,
    tableRows,
    tableColumnNamesEnrichedGenes,
    isPastBreakpoint,
  ]);

  const tableUnavailableComponent = (
    <TableUnavailableContainer>
      <TableUnavailableHeader>
        No {activeTable ? "computational" : "canonical"} marker genes
      </TableUnavailableHeader>
      <TableUnavailableDescription>
        {activeTable
          ? getEmptyComputationalMarkerGenesTableUIMessageDetail(
              allFilteredByLowMarkerScore
            )
          : getEmptyCanonicalMarkerGenesTableUIMessageDetail()}
      </TableUnavailableDescription>
    </TableUnavailableContainer>
  );

  const handlePageChange = (
    _event: React.ChangeEvent<unknown>,
    page: number
  ) => {
    setPage(page);
    track(EVENTS.CG_MARKER_GENE_PAGINATION_CLICKED, {
      page: page,
      type: activeTable ? "computational" : "canonical",
      cell_type: cellTypeName,
      inSideBar: !!cellInfoCellType,
    });
  };

  return (
    <MarkerGeneTableWrapper>
      <TableTitleWrapper>
        <TableTitleOuterWrapper>
          <TableTitleInnerWrapper columnGap={4}>
            <TableTitle>Marker Genes</TableTitle>
          </TableTitleInnerWrapper>
          <TableTitleInnerWrapper>
            {!skinnyMode && !cellInfoCellType && (
              <Link
                url={`${ROUTES.WHERE_IS_MY_GENE}?genes=${genesForShareUrl}&ver=2`}
                label="Open in Gene Expression"
                onClick={() => {
                  track(EVENTS.CG_OPEN_IN_WMG_CLICKED, {
                    type: activeTable ? "computational" : "canonical",
                  });
                }}
              />
            )}
          </TableTitleInnerWrapper>
        </TableTitleOuterWrapper>
      </TableTitleWrapper>
      <TableTitleOuterWrapper>
        <TableSelectorRow>
          <FlexRow>
            <TableSelectorButton
              data-testid={CELL_GUIDE_CARD_ENRICHED_GENES_TABLE_SELECTOR}
              isActive={activeTable === 1}
              onClick={() => {
                setPage(1);
                setActiveTable(1);
                track(EVENTS.CG_COMPUTATIONAL_TAB_CLICKED, {
                  cell_type: cellTypeName,
                  inSideBar: !!cellInfoCellType,
                });
              }}
            >
              Computational
            </TableSelectorButton>
          </FlexRow>
          <FlexRow>
            <TableSelectorButton
              data-testid={
                CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE_SELECTOR
              }
              isActive={activeTable === 0}
              onClick={() => {
                setPage(1);
                setActiveTable(0);
                track(EVENTS.CG_CANONICAL_TAB_CLICKED, {
                  cell_type: cellTypeName,
                  inSideBar: !!cellInfoCellType,
                });
              }}
            >
              Canonical
            </TableSelectorButton>
          </FlexRow>
        </TableSelectorRow>
      </TableTitleOuterWrapper>
      {tableRows.length > 0 ? (
        <div ref={containerRef}>
          {tableComponent}
          <MarkerGenePagination>
            <Pagination
              count={pageCount}
              page={page}
              onChange={handlePageChange}
            />
            {activeTable === 0 ? (
              <MarkerGeneInfo>
                Source: HuBMAP
                <HelpTooltip
                  skinnyMode={skinnyMode}
                  title="Canonical Marker Genes"
                  setTooltipContent={setTooltipContent}
                  dark
                  placement="top-end"
                  buttonDataTestId={MARKER_GENES_CANONICAL_TOOLTIP_TEST_ID}
                  text={canonicalMarkerGenesTooltipComponent}
                />
              </MarkerGeneInfo>
            ) : (
              <MarkerGeneInfo>
                Source: CZI
                <HelpTooltip
                  skinnyMode={skinnyMode}
                  title="Computational Marker Genes"
                  setTooltipContent={setTooltipContent}
                  dark
                  placement="top-end"
                  buttonDataTestId={MARKER_GENES_COMPUTATIONAL_TOOLTIP_TEST_ID}
                  text={enrichedGenesTooltipComponent}
                />
              </MarkerGeneInfo>
            )}
          </MarkerGenePagination>
        </div>
      ) : (
        tableUnavailableComponent
      )}
    </MarkerGeneTableWrapper>
  );
};

export default MarkerGeneTables;
