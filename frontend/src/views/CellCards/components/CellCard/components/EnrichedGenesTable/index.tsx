import React, { useEffect, useMemo, useState } from "react";
import {
  TableTitle,
  TableTitleWrapper,
  TableUnavailableContainer,
  TableUnavailableHeader,
  TableUnavailableDescription,
} from "../common/style";
import Link from "../common/Link";
import { TableTitleInnerWrapper, TableTitleOuterWrapper } from "./style";
import Table from "../Table";
import { HIGHLY_EXPRESSED_GENES_SECTION_ID } from "../CellCardSidebar";
import DropdownSelect from "../common/DropdownSelect";
import { SelectChangeEvent } from "@mui/material/Select";
import { useEnrichedGenes } from "src/common/queries/cellCards";
import HelpTooltip from "../common/HelpTooltip";
import { ROUTES } from "src/common/constants/routes";

interface TableRow {
  symbol: string;
  name: string;
  me: string;
  pc: string;
}
const tableColumns: Array<keyof TableRow> = ["symbol", "name", "me", "pc"];

const tableColumnNames: Record<keyof TableRow, string> = {
  symbol: "Symbol",
  name: "Name",
  me: "Expression Score",
  pc: "% of Cells",
};

interface Props {
  cellTypeId: string;
}

const EnrichedGenesTable = ({ cellTypeId }: Props) => {
  const { data: genes } = useEnrichedGenes(cellTypeId);
  const [selectedOrganism, setSelectedOrganism] = useState("");

  const uniqueOrganisms = useMemo(() => {
    if (!genes) return [];
    const organisms = new Set<string>();
    for (const markerGene of genes) {
      organisms.add(markerGene.organism);
    }
    if (!selectedOrganism && Array.from(organisms).includes("Homo sapiens"))
      setSelectedOrganism("Homo sapiens");
    else if (!selectedOrganism)
      setSelectedOrganism(Array.from(organisms).at(0) ?? "");

    const organismsArray = Array.from(organisms).sort((a, b) => {
      if (a === "Homo sapiens") return -1;
      if (b === "Homo sapiens") return 1;
      return a.localeCompare(b);
    });
    return organismsArray;
  }, [genes, cellTypeId]);

  const tableRows: TableRow[] = useMemo(() => {
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

  useEffect(() => {
    return () => {
      setSelectedOrganism("");
    };
  }, []);

  // const genesForShareUrl = tableRows.map((row) => row.symbol).join("%2C");

  const handleChange = (event: SelectChangeEvent) => {
    setSelectedOrganism(event.target.value as string);
  };

  return (
    <>
      <TableTitleWrapper id={HIGHLY_EXPRESSED_GENES_SECTION_ID}>
        <TableTitleOuterWrapper>
          <TableTitleInnerWrapper columnGap={4}>
            <TableTitle>Highly Expressed Genes</TableTitle>
            <HelpTooltip
              text={
                <>
                  {
                    "The marker genes listed below are computationally derived from the "
                  }
                  <Link label={"CELLxGENE corpus"} url={ROUTES.DATASETS} />
                  {
                    ". They are computed utilizing the same methodology as featured in our "
                  }
                  <Link
                    label={
                      "Find Marker Genes feature from the Gene Expression application"
                    }
                    url={ROUTES.FMG_DOCS}
                  />
                  {"."}
                </>
              }
            />
          </TableTitleInnerWrapper>
          {tableRows.length > 0 && (
            <TableTitleInnerWrapper>
              <DropdownSelect
                handleChange={handleChange}
                options={uniqueOrganisms}
                selectedOption={selectedOrganism}
              />
              <Link
                url={`https://cellxgene.cziscience.com/gene-expression`}
                label="Open in Gene Expression"
              />
            </TableTitleInnerWrapper>
          )}
        </TableTitleOuterWrapper>
      </TableTitleWrapper>
      {tableRows.length > 0 ? (
        <Table<TableRow>
          columns={tableColumns}
          rows={tableRows}
          columnIdToName={tableColumnNames}
        />
      ) : (
        <TableUnavailableContainer>
          <TableUnavailableHeader>
            No highly expressed genes
          </TableUnavailableHeader>
          <TableUnavailableDescription>
            Highly expressed genes for this cell type are unavailable at this
            time.
          </TableUnavailableDescription>
        </TableUnavailableContainer>
      )}
    </>
  );
};
export default EnrichedGenesTable;
