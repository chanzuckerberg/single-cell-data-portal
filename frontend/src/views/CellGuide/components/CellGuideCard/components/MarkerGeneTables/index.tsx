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
} from "./style";
import Table from "../common/Table";
import DropdownSelect from "../common/DropdownSelect";
import { SelectChangeEvent } from "@mui/material/Select";
import { Pagination } from "@mui/material";
import {
  useCanonicalMarkers,
  useEnrichedGenes,
} from "src/common/queries/cellGuide";
import HelpTooltip from "../common/HelpTooltip";
import { ROUTES } from "src/common/constants/routes";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { FMG_GENE_STRENGTH_THRESHOLD } from "src/views/WheresMyGene/common/constants";
import { CENSUS_LINK } from "src/components/Header/components/Nav";

export const CELL_GUIDE_CARD_MARKER_GENES_TABLE_DROPDOWN_ORGANISM =
  "cell-guide-card-marker-genes-table-dropdown-organism";
export const CELL_GUIDE_CARD_MARKER_GENES_TABLE_DROPDOWN_ORGAN =
  "cell-guide-card-marker-genes-table-dropdown-organ";
export const CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE =
  "cell-guide-card-canonical-marker-genes-table";
export const CELL_GUIDE_CARD_ENRICHED_GENES_TABLE =
  "cell-guide-card-enriched-genes-table";
export const CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE_SELECTOR =
  "cell-guide-card-canonical-marker-genes-table-selector";
export const CELL_GUIDE_CARD_ENRICHED_GENES_TABLE_SELECTOR =
  "cell-guide-card-enriched-genes-table-selector";

export const EXPRESSION_SCORE_TOOLTIP_TEST_ID =
  "cell-guide-card-expression-score-tooltip";
export const PERCENT_OF_CELLS_TOOLTIP_TEST_ID =
  "cell-guide-card-percent-of-cells-tooltip";
export const MARKER_SCORE_TOOLTIP_TEST_ID =
  "cell-guide-card-marker-score-tooltip";

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
        <HelpTooltipWrapper
          buttonDataTestId={MARKER_SCORE_TOOLTIP_TEST_ID}
          content={
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
      <HelpTooltipWrapper
        buttonDataTestId={EXPRESSION_SCORE_TOOLTIP_TEST_ID}
        content={
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
      <HelpTooltipWrapper
        buttonDataTestId={PERCENT_OF_CELLS_TOOLTIP_TEST_ID}
        content={
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
const tableColumnNamesCanonicalGenes: Record<
  keyof TableRowCanonicalGenes,
  string
> = {
  symbol: "Symbol",
  name: "Name",
  references: "References",
};

type TableRow = (TableRowEnrichedGenes | TableRowCanonicalGenes) & {
  symbolId: string;
};

interface Props {
  cellTypeId: string;
  setGeneInfoGene: React.Dispatch<React.SetStateAction<string | null>>;
  cellTypeName: string;
}

const ROWS_PER_PAGE = 10;

export const MARKER_GENES_COMPUTATIONAL_TOOLTIP_TEST_ID =
  "marker-genes-computational-help-tooltip";
export const MARKER_GENES_CANONICAL_TOOLTIP_TEST_ID =
  "marker-genes-canonical-help-tooltip";

const MarkerGeneTables = ({
  cellTypeId,
  cellTypeName,
  setGeneInfoGene,
}: Props) => {
  const [selectedOrganism, setSelectedOrganism] = useState("");
  // 0 is canonical marker genes, 1 is computational marker genes
  const [activeTable, setActiveTable] = useState(0);
  const [selectedOrgan, setSelectedOrgan] = useState("");
  const [page, setPage] = useState(1);

  let uniqueOrganisms = ["Homo sapiens"];
  let uniqueOrgans: string[] = ["All Tissues"];
  let tableRows: TableRow[];
  if (activeTable) {
    const { data: genes } = useEnrichedGenes(cellTypeId);
    uniqueOrganisms = useMemo(() => {
      if (!genes) return [];
      const organisms = new Set<string>();
      for (const markerGene of genes) {
        organisms.add(markerGene.organism);
      }
      if (!selectedOrganism && Array.from(organisms).includes("Homo sapiens"))
        setSelectedOrganism("Homo sapiens");
      else if (!selectedOrganism)
        setSelectedOrganism(Array.from(organisms).at(0) ?? "");
      // All Tissues always selected
      setSelectedOrgan("All Tissues");

      return Array.from(organisms).sort((a, b) => {
        if (a === "Homo sapiens") return -1;
        if (b === "Homo sapiens") return 1;
        return a.localeCompare(b);
      });
    }, [genes, cellTypeId]);

    tableRows = useMemo(() => {
      if (!genes) return [];
      const rows: TableRow[] = [];
      for (const markerGene of genes) {
        const { pc, me, name, symbol, organism, marker_score } = markerGene;
        if (organism !== selectedOrganism) continue;
        if (marker_score < FMG_GENE_STRENGTH_THRESHOLD) continue;
        rows.push({
          symbolId: symbol,
          symbol: (
            <>
              {symbol}{" "}
              <ButtonIcon
                aria-label={`display gene info for ${symbol}`}
                sdsIcon="infoCircle"
                sdsSize="small"
                sdsType="secondary"
                onClick={() => setGeneInfoGene(symbol.toUpperCase())}
              />
            </>
          ),
          name,
          marker_score: marker_score.toFixed(2),
          me: me.toFixed(2),
          pc: (pc * 100).toFixed(1),
        });
      }
      return rows;
    }, [genes, selectedOrganism, cellTypeId]);
  } else {
    const { data: genes } = useCanonicalMarkers(cellTypeId);
    uniqueOrgans = useMemo(() => {
      if (!genes) return [];
      const organs = new Set<string>();
      for (const markerGene of genes) {
        organs.add(markerGene.tissue);
      }
      // All Tissues selecteed by default
      if (!selectedOrgan) setSelectedOrgan("All Tissues");
      // Homo sapiens always selected
      setSelectedOrganism("Homo sapiens");
      const uniqueOrgans = Array.from(organs);

      // put "All Tissues" first in the array
      return uniqueOrgans.sort((a, b) => {
        if (a === "All Tissues") return -1;
        if (b === "All Tissues") return 1;
        return a.localeCompare(b);
      });
    }, [genes, cellTypeId]);

    tableRows = useMemo(() => {
      if (!genes) return [];
      const rows: (TableRow & {
        numReferences: number;
      })[] = [];

      const publicationTitlesToIndex = new Map();
      let index = 0;
      for (const markerGene of genes) {
        if (markerGene.tissue !== selectedOrgan) continue;
        const publicationTitles = markerGene.publication_titles.split(";;");
        for (let i = 0; i < publicationTitles.length; i += 1) {
          if (
            !publicationTitlesToIndex.has(publicationTitles[i]) &&
            publicationTitles[i] !== ""
          ) {
            publicationTitlesToIndex.set(publicationTitles[i], index);
            index += 1;
          }
        }
      }
      for (const markerGene of genes) {
        const { tissue, publication, publication_titles, symbol, name } =
          markerGene;

        if (tissue !== selectedOrgan) continue;

        // multiple publications for a single gene are joined by ";;"
        let publications = Array.from(new Set(publication.split(";;")));
        let publicationTitles = Array.from(
          new Set(publication_titles.split(";;"))
        );

        const sortedPublicationsAndTitles = publications
          .map((pub, i) => [pub, publicationTitles[i]])
          .sort((a, b) => {
            return (
              publicationTitlesToIndex.get(a[1]) -
              publicationTitlesToIndex.get(b[1])
            );
          })
          .filter((publicationTitle, index) => {
            return publicationTitle && publications[index];
          });

        publications = sortedPublicationsAndTitles.map((pub) => pub[0]);
        publicationTitles = sortedPublicationsAndTitles.map((pub) => pub[1]);

        const publicationLinks = (
          <PublicationLinkWrapper>
            {publicationTitles.map((publicationTitle, index) => {
              if (publicationTitle && publications[index]) {
                const referenceIndexLabel =
                  publicationTitlesToIndex.get(publicationTitle) + 1;
                return (
                  <Tooltip
                    key={`${publications[index]}-${index}-tooltip`}
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
                    <span key={`${publications[index]}-${index}-span`}>
                      <Link
                        key={`${publications[index]}-${index}`}
                        label={`[${referenceIndexLabel}]`}
                        url={`https://doi.org/${publications[index]}`}
                      />
                    </span>
                  </Tooltip>
                );
              }
            })}
          </PublicationLinkWrapper>
        );

        rows.push({
          name,
          symbolId: symbol,
          symbol: (
            <>
              {symbol}{" "}
              <ButtonIcon
                aria-label={`display gene info for ${symbol}`}
                sdsIcon="infoCircle"
                sdsSize="small"
                sdsType="secondary"
                onClick={() => setGeneInfoGene(symbol.toUpperCase())}
              />
            </>
          ),
          references: publicationLinks,
          numReferences: sortedPublicationsAndTitles.length,
        });
      }

      // Sort rows by number of references
      if (rows.length) {
        rows.sort((a, b) => {
          return b.numReferences - a.numReferences;
        });
      }

      return rows;
    }, [genes, selectedOrgan, setGeneInfoGene]);
  }

  useEffect(() => {
    return () => {
      // (alec) when the component unmounts, reset the organism/organ to its initial state.
      // not all cell types may have homo sapiens as a valid option so we need to invoke the conditional
      // logic above for setting the initial organism/organ.
      setSelectedOrganism("");
      setSelectedOrgan("");
      setActiveTable(0);
      setPage(1);
    };
  }, []);

  // Handle cell type change, set marker genes table page back to 1
  useEffect(() => {
    setPage(1);
  }, [cellTypeId]);

  const genesForShareUrl = `${tableRows
    .map((row) => row.symbolId)
    .join("%2C")}&cellTypes=${cellTypeName.replace(" ", "+")}`;

  const handleChangeOrganism = (event: SelectChangeEvent<unknown>) => {
    setSelectedOrganism(event.target.value as string);
  };
  const handleChangeOrgan = (event: SelectChangeEvent<unknown>) => {
    setSelectedOrgan(event.target.value as string);
  };
  const enrichedGenesTooltipComponent = (
    <div>
      {"Computational marker genes are derived from the "}
      <Link label={"CELLxGENE Census"} url={CENSUS_LINK} />
      {". They are computed utilizing the same methodology as featured in our "}
      <Link
        label={"Find Marker Genes feature from the Gene Expression application"}
        url={ROUTES.FMG_DOCS}
      />
      {"."}
    </div>
  );
  const canonicalMarkerGenesTooltipComponent = (
    <div>
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
      <br />
      <br />
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
    </div>
  );

  const pageCount = Math.ceil(tableRows.length / ROWS_PER_PAGE);
  const tableComponent = activeTable ? (
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
  const tableUnavailableComponent = (
    <TableUnavailableContainer>
      <TableUnavailableHeader>
        No {activeTable ? "computational" : "canonical"} marker genes
      </TableUnavailableHeader>
      <TableUnavailableDescription>
        {activeTable ? "Computational" : "Canonical"} marker genes for this cell
        type are unavailable at this time.
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
          {tableRows.length > 0 && (
            <TableTitleInnerWrapper>
              {activeTable === 1 && (
                <DropdownSelect
                  handleChange={handleChangeOrganism}
                  options={uniqueOrganisms}
                  selectedOption={selectedOrganism}
                  testId={CELL_GUIDE_CARD_MARKER_GENES_TABLE_DROPDOWN_ORGANISM}
                />
              )}
              {activeTable === 0 && (
                <DropdownSelect
                  handleChange={handleChangeOrgan}
                  options={uniqueOrgans}
                  selectedOption={selectedOrgan}
                  testId={CELL_GUIDE_CARD_MARKER_GENES_TABLE_DROPDOWN_ORGAN}
                />
              )}
              <Link
                url={`${ROUTES.WHERE_IS_MY_GENE}?genes=${genesForShareUrl}&ver=2`}
                label="Open in Gene Expression"
              />
            </TableTitleInnerWrapper>
          )}
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
              Canonical (HuBMAP)
            </TableSelectorButton>
            <HelpTooltip
              buttonDataTestId={MARKER_GENES_CANONICAL_TOOLTIP_TEST_ID}
              text={canonicalMarkerGenesTooltipComponent}
            />
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
              Computational (CZI)
            </TableSelectorButton>
            <HelpTooltip
              buttonDataTestId={MARKER_GENES_COMPUTATIONAL_TOOLTIP_TEST_ID}
              text={enrichedGenesTooltipComponent}
            />
          </FlexRow>
        </TableSelectorRow>
      </TableTitleOuterWrapper>
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

interface HelpTooltipWrapperProps {
  buttonDataTestId: string;
  content: ReactElement;
}
function HelpTooltipWrapper({
  buttonDataTestId,
  content,
}: HelpTooltipWrapperProps) {
  return (
    <HelpTooltip dark buttonDataTestId={buttonDataTestId} text={content} />
  );
}

export default MarkerGeneTables;
