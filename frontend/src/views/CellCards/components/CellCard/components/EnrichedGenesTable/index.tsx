import React, { ReactElement, useMemo, useState } from "react";
import {
  TableTitle,
  TableTitleWrapper,
  WmgLink,
  TableUnavailableContainer,
  TableUnavailableHeader,
  TableUnavailableDescription,
} from "../common/style";
import { TableTitleInnerWrapper, TableTitleOuterWrapper } from "./style";
import Table from "../Table";
import { HIGHLY_EXPRESSED_GENES_SECTION_ID } from "../CellCardSidebar";
import Link from "../common/Link";
import DropdownSelect from "../common/DropdownSelect";
import { SelectChangeEvent } from "@mui/material/Select";
import { useEnrichedGenes } from "src/common/queries/cellCards";

interface TableRow {
  symbol: ReactElement;
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
      if (!selectedOrganism) setSelectedOrganism(markerGene.organism);
    }
    return Array.from(organisms);
  }, [genes]);

  const tableRows: TableRow[] = useMemo(() => {
    if (!genes) return [];
    const rows = [];
    for (const markerGene of genes) {
      const { pc, me, name, symbol, organism } = markerGene;
      if (organism !== selectedOrganism) continue;
      rows.push({
        symbol: (
          <Link
            title={`${symbol}`}
            url={`https://www.genecards.org/cgi-bin/carddisp.pl?gene=${markerGene.symbol}`}
          />
        ),
        name,
        me: me.toFixed(2),
        pc: (pc * 100).toFixed(1),
      });
    }
    return rows;
  }, [genes, selectedOrganism]);

  const genesForShareUrl = tableRows.map((row) => row.symbol).join("%2C");

  const handleChange = (event: SelectChangeEvent) => {
    setSelectedOrganism(event.target.value as string);
  };

  return (
    <>
      <TableTitleWrapper id={HIGHLY_EXPRESSED_GENES_SECTION_ID}>
        <TableTitleOuterWrapper>
          <TableTitleInnerWrapper>
            <TableTitle>Highly Expressed Genes</TableTitle>
            {tableRows.length > 0 && (
              <DropdownSelect
                handleChange={handleChange}
                options={uniqueOrganisms}
                selectedOption={selectedOrganism}
              />
            )}
          </TableTitleInnerWrapper>
          {tableRows.length > 0 && (
            <WmgLink
              href={`https://cellxgene.cziscience.com/gene-expression?tissues=lung&genes=${genesForShareUrl}&ver=2`}
              target="_blank"
            >
              Open in Gene Expression
            </WmgLink>
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
