import React, { ReactElement, useEffect, useMemo, useState } from "react";
import { Tooltip } from "@czi-sds/components";
import {
  TableTitle,
  TableTitleWrapper,
  TableUnavailableContainer,
  TableUnavailableHeader,
  TableUnavailableDescription,
  TableTitleInnerWrapper,
  StyledDivider,
} from "../common/style";
import Link from "../common/Link";
import {
  PublicationLinkWrapper,
  TableSelectorButton,
  TableSelectorRow,
  TableTitleOuterWrapper,
} from "./style";
import Table from "../common/Table";
import DropdownSelect from "../common/DropdownSelect";
import { SelectChangeEvent } from "@mui/material/Select";
import {
  useCanonicalMarkers,
  useEnrichedGenes,
} from "src/common/queries/cellCards";
import HelpTooltip from "../common/HelpTooltip";
import { ROUTES } from "src/common/constants/routes";

export const CELL_CARD_MARKER_GENES_TABLE_DROPDOWN_ORGANISM =
  "cell-card-marker-genes-table-dropdown-organism";
export const CELL_CARD_MARKER_GENES_TABLE_DROPDOWN_ORGAN =
  "cell-card-marker-genes-table-dropdown-organ";
export const CELL_CARD_CANONICAL_MARKER_GENES_TABLE =
  "cell-card-canonical-marker-genes-table";
export const CELL_CARD_ENRICHED_GENES_TABLE = "cell-card-enriched-genes-table";
export const CELL_CARD_CANONICAL_MARKER_GENES_TABLE_SELECTOR =
  "cell-card-canonical-marker-genes-table-selector";
export const CELL_CARD_ENRICHED_GENES_TABLE_SELECTOR =
  "cell-card-enriched-genes-table-selector";

interface TableRowEnrichedGenes {
  symbol: string;
  name: string;
  me: string;
  pc: string;
}
const tableColumnsEnrichedGenes: Array<keyof TableRowEnrichedGenes> = [
  "symbol",
  "name",
  "me",
  "pc",
];

const tableColumnNamesEnrichedGenes: Record<
  keyof TableRowEnrichedGenes,
  string
> = {
  symbol: "Symbol",
  name: "Name",
  me: "Expression Score",
  pc: "% of Cells",
};

interface TableRowCanonicalGenes {
  symbol: string;
  name: string;
  references: ReactElement | string;
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

type TableRow = TableRowEnrichedGenes | TableRowCanonicalGenes;

interface Props {
  cellTypeId: string;
}

const MarkerGeneTables = ({ cellTypeId }: Props) => {
  const [selectedOrganism, setSelectedOrganism] = useState("");
  const [activeTable, setActiveTable] = useState(0);
  const [selectedOrgan, setSelectedOrgan] = useState("");

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

      const organismsArray = Array.from(organisms).sort((a, b) => {
        if (a === "Homo sapiens") return -1;
        if (b === "Homo sapiens") return 1;
        return a.localeCompare(b);
      });
      return organismsArray;
    }, [genes, cellTypeId]);

    tableRows = useMemo(() => {
      if (!genes) return [];
      const rows = [];
      for (const markerGene of genes) {
        const { pc, me, name, symbol, organism } = markerGene;
        if (organism !== selectedOrganism) continue;
        rows.push({
          symbol: symbol,
          name,
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
        organs.add(markerGene.tissue_general);
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
      const rows = [];

      for (const markerGene of genes) {
        if (markerGene.tissue_general !== selectedOrgan) continue;
        // multiple publications for a single gene are joined by ";;"
        const publications = markerGene.publication.split(";;");
        const publicationTitles = markerGene.publication_titles.split(";;");
        const publicationLinks = (
          <PublicationLinkWrapper>
            {publications.map((publication, index) => {
              if (publication && publicationTitles[index]) {
                return (
                  <Tooltip
                    key={`${publication}-${index}-tooltip`}
                    placement="right"
                    width="default"
                    arrow={false}
                    title={publicationTitles[index]}
                    leaveDelay={0}
                  >
                    <span key={`${publication}-${index}-span`}>
                      <Link
                        key={`${publication}-${index}`}
                        label={`[${index + 1}]`}
                        url={`https://doi.org/${publication}`}
                      />
                    </span>
                  </Tooltip>
                );
              }
            })}
          </PublicationLinkWrapper>
        );

        rows.push({
          symbol: markerGene.symbol,
          name: markerGene.name,
          references: publicationLinks,
        });
      }
      return rows;
    }, [genes, selectedOrgan, cellTypeId]);
  }

  useEffect(() => {
    return () => {
      // (alec) when the component unmounts, reset the organism/organ to its initial state.
      // not all cell types may have homo sapiens as a valid option so we need to invoke the conditional
      // logic above for setting the initial organism/organ.
      setSelectedOrganism("");
      setSelectedOrgan("");
    };
  }, []);

  const genesForShareUrl = tableRows.map((row) => row.symbol).join("%2C");

  const handleChangeOrganism = (event: SelectChangeEvent<unknown>) => {
    setSelectedOrganism(event.target.value as string);
  };
  const handleChangeOrgan = (event: SelectChangeEvent<unknown>) => {
    setSelectedOrgan(event.target.value as string);
  };
  const enrichedGenesTooltipComponent = (
    <div>
      {"The marker genes listed below are computationally derived from the "}
      <Link label={"CELLxGENE corpus"} url={ROUTES.DATASETS} />
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
        "The below marker genes and associated publications were derived from the "
      }
      <Link
        label={"Anatomical Structures, Cell Types and Biomarkers (ASCT+B)"}
        url={"https://humanatlas.io/asctb-tables"}
      />
      {
        " tables. The tables are authored and reviewed by an international team of anatomists, pathologists, physicians, and other experts."
      }
      <br />
      <br />
      <i>
        Quardokus, Ellen, Bruce W. Herr II, Lisel Record, Katy Börner. 2022.
        HuBMAP ASCT+B Tables. Accessed May 16, 2023.
      </i>
    </div>
  );
  const tableComponent = activeTable ? (
    <Table<TableRowEnrichedGenes>
      columns={tableColumnsEnrichedGenes}
      rows={tableRows as TableRowEnrichedGenes[]}
      columnIdToName={tableColumnNamesEnrichedGenes}
    />
  ) : (
    <Table<TableRowCanonicalGenes>
      columns={tableColumnsCanonicalGenes}
      rows={tableRows as TableRowCanonicalGenes[]}
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

  return (
    <div
      data-testid={
        activeTable
          ? CELL_CARD_ENRICHED_GENES_TABLE
          : CELL_CARD_CANONICAL_MARKER_GENES_TABLE
      }
    >
      <TableTitleWrapper>
        <TableTitleOuterWrapper>
          <TableTitleInnerWrapper columnGap={4}>
            <TableTitle>Marker Genes</TableTitle>
            <HelpTooltip
              text={
                activeTable
                  ? enrichedGenesTooltipComponent
                  : canonicalMarkerGenesTooltipComponent
              }
            />
          </TableTitleInnerWrapper>
          {tableRows.length > 0 && (
            <TableTitleInnerWrapper>
              <DropdownSelect
                handleChange={handleChangeOrganism}
                options={uniqueOrganisms}
                selectedOption={selectedOrganism}
                testId={CELL_CARD_MARKER_GENES_TABLE_DROPDOWN_ORGANISM}
              />
              <DropdownSelect
                handleChange={handleChangeOrgan}
                options={uniqueOrgans}
                selectedOption={selectedOrgan}
                testId={CELL_CARD_MARKER_GENES_TABLE_DROPDOWN_ORGAN}
              />
              <Link
                url={`${ROUTES.WHERE_IS_MY_GENE}?genes=${genesForShareUrl}&ver=2`}
                label="Open in Gene Expression"
              />
            </TableTitleInnerWrapper>
          )}
        </TableTitleOuterWrapper>
      </TableTitleWrapper>
      <TableSelectorRow>
        <TableSelectorButton
          data-testid={CELL_CARD_CANONICAL_MARKER_GENES_TABLE_SELECTOR}
          isActive={activeTable === 0}
          onClick={() => setActiveTable(0)}
        >
          Canonical
        </TableSelectorButton>
        <TableSelectorButton
          data-testid={CELL_CARD_ENRICHED_GENES_TABLE_SELECTOR}
          isActive={activeTable === 1}
          onClick={() => setActiveTable(1)}
        >
          Computational
        </TableSelectorButton>
      </TableSelectorRow>
      <StyledDivider />
      {tableRows.length > 0 ? tableComponent : tableUnavailableComponent}
    </div>
  );
};
export default MarkerGeneTables;
