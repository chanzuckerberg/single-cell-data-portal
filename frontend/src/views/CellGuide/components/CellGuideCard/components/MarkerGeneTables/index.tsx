/* eslint-disable react/no-unescaped-entities */
import React, {
  ReactElement,
  ReactNode,
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
import { useComputationalMarkerGenesTableRowsAndFilters } from "./hooks/computational_markers";
import { useCanonicalMarkerGenesTableRowsAndFilters } from "./hooks/canonical_markers";
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
} from "src/views/CellGuide/components/CellGuideCard/components/MarkerGeneTables/constants";
import { FMG_GENE_STRENGTH_THRESHOLD } from "src/views/WheresMyGene/common/constants";

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
  name: string;
  marker_score: string;
  me: string;
  pc: string;
}
const tableColumnsEnrichedGenes: Array<keyof TableRowEnrichedGenes> = [
  "symbol",
  "name",
  "marker_score",
  "me",
  "pc",
];

// Computational marker gene table column names
const tableColumnNamesEnrichedGenes: Record<
  keyof TableRowEnrichedGenes,
  ReactElement | string
> = {
  symbol: "Symbol",
  name: "Name",
  marker_score: (
    <div>
      <StyledHeadCellContent>
        Marker Score
        <HelpTooltip
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
                Marker genes are highly and uniquely expressed in the cell type
                relative to all other cell types.
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
      Expression Score
      <HelpTooltip
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
};

// Canonical marker gene table types
interface TableRowCanonicalGenes {
  symbol: ReactNode;
  name: string;
  references: ReactNode;
}
const tableColumnsCanonicalGenes: Array<keyof TableRowCanonicalGenes> = [
  "symbol",
  "name",
  "references",
];

// Canonical marker gene table column names
const tableColumnNamesCanonicalGenes: Record<
  keyof TableRowCanonicalGenes,
  string
> = {
  symbol: "Symbol",
  name: "Name",
  references: "References",
};

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
  organName: string;
  organId: string;
  organismName: string;
}

const MarkerGeneTables = ({
  cellTypeId,
  cellTypeName,
  setGeneInfoGene,
  organName,
  organId,
  organismName,
}: Props) => {
  // 0 is canonical marker genes, 1 is computational marker genes
  const [activeTable, setActiveTable] = useState(0);
  const [computationalMarkerGenes, setComputationalMarkerGenes] =
    useState<ComputationalMarkersQueryResponse>([]);
  const [canonicalMarkerGenes, setCanonicalMarkerGenes] =
    useState<CanonicalMarkersQueryResponse>([]);

  const { data: enrichedGenes } = useComputationalMarkers(cellTypeId);
  const { data: canonicalMarkers } = useCanonicalMarkers(cellTypeId);
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
      selectedOrganLabel: organName,
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

  const tableRows: TableRow[] = useMemo(
    () =>
      activeTable
        ? computationalMarkerGeneTableData.map((row) => ({
            ...row,
            symbolId: row.symbol,
            symbol: (
              <>
                {row.symbol}{" "}
                <ButtonIcon
                  aria-label={`display gene info for ${row.symbol}`}
                  sdsIcon="infoCircle"
                  sdsSize="small"
                  sdsType="secondary"
                  onClick={() => setGeneInfoGene(row.symbol.toUpperCase())}
                />
              </>
            ),
          }))
        : canonicalMarkerGeneTableData.map((row) => ({
            ...row,
            symbolId: row.symbol,
            symbol: (
              <>
                {row.symbol}{" "}
                <ButtonIcon
                  aria-label={`display gene info for ${row.symbol}`}
                  sdsIcon="infoCircle"
                  sdsSize="small"
                  sdsType="secondary"
                  onClick={() => setGeneInfoGene(row.symbol.toUpperCase())}
                />
              </>
            ),
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
                      return (
                        <Tooltip
                          key={`${row.referenceData.publications[index]}-${index}-tooltip`}
                          placement="top"
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
                            key={`${row.referenceData.publications[index]}-${index}-span`}
                          >
                            <Link
                              key={`${row.referenceData.publications[index]}-${index}`}
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
          })),
    [
      activeTable,
      canonicalMarkerGeneTableData,
      computationalMarkerGeneTableData,
      setGeneInfoGene,
    ]
  );

  const genesForShareUrl = `${tableRows
    .map((row) => row.symbolId)
    .join("%2C")}&cellTypes=${cellTypeName.replace(" ", "+")}`;

  const pageCount = Math.ceil(tableRows.length / ROWS_PER_PAGE);
  const tableComponent = useMemo(
    () =>
      activeTable ? (
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
      ),
    [activeTable, page, tableRows]
  );

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
  };

  return (
    <div>
      <TableTitleWrapper>
        <TableTitleOuterWrapper>
          <TableTitleInnerWrapper columnGap={4}>
            <TableTitle>Marker Genes</TableTitle>
          </TableTitleInnerWrapper>
          <TableTitleInnerWrapper>
            <Link
              url={`${ROUTES.WHERE_IS_MY_GENE}?genes=${genesForShareUrl}&ver=2`}
              label="Open in Gene Expression"
            />
          </TableTitleInnerWrapper>
        </TableTitleOuterWrapper>
      </TableTitleWrapper>
      <TableTitleOuterWrapper>
        <TableSelectorRow>
          <FlexRow>
            <TableSelectorButton
              data-testid={
                CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE_SELECTOR
              }
              isActive={activeTable === 0}
              onClick={() => {
                setPage(1);
                setActiveTable(0);
                track(EVENTS.CG_CANONICAL_TAB_CLICKED);
              }}
            >
              Canonical
            </TableSelectorButton>
          </FlexRow>
          <FlexRow>
            <TableSelectorButton
              data-testid={CELL_GUIDE_CARD_ENRICHED_GENES_TABLE_SELECTOR}
              isActive={activeTable === 1}
              onClick={() => {
                setPage(1);
                setActiveTable(1);
                track(EVENTS.CG_COMPUTATIONAL_TAB_CLICKED);
              }}
            >
              Computational
            </TableSelectorButton>
          </FlexRow>
        </TableSelectorRow>
      </TableTitleOuterWrapper>
      {tableRows.length > 0 ? (
        <div>
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
    </div>
  );
};

export default MarkerGeneTables;
