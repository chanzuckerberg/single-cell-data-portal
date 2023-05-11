import React, { ReactElement } from "react";
import { TableTitle, TableTitleWrapper } from "../common/style";
import { WmgLink } from "../CanonicalMarkerGeneTable/style";
import Table from "../Table";
import { allEnrichedGenes } from "src/views/CellCards/common/fixtures";
import { HIGHLY_EXPRESSED_GENES_SECTION_ID } from "../../../CellCardSidebar";

const Link = ({ title, url }: { title: string; url: string }) => {
  return (
    <a href={url} target="_blank">
      {title}
    </a>
  );
};

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
  const tableRows: TableRow[] = [];
  if (cellTypeId in allEnrichedGenes) {
    const genes = allEnrichedGenes[cellTypeId as keyof typeof allEnrichedGenes];
    for (const markerGene of genes) {
      const { pc, me, name, symbol } = markerGene;

      tableRows.push({
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
  }

  const genesForShareUrl = tableRows.map((row) => row.symbol).join("%2C");

  return (
    <>
      <TableTitleWrapper id={HIGHLY_EXPRESSED_GENES_SECTION_ID}>
        <TableTitle>Highly Expressed Genes</TableTitle>
        {tableRows.length > 0 && (
          <WmgLink
            href={`https://cellxgene.cziscience.com/gene-expression?tissues=lung&genes=${genesForShareUrl}&ver=2`}
            target="_blank"
          >
            Open in Gene Expression
          </WmgLink>
        )}
      </TableTitleWrapper>
      {tableRows.length && (
        <Table<TableRow>
          columns={tableColumns}
          rows={tableRows}
          columnIdToName={tableColumnNames}
        />
      )}
    </>
  );
};
export default EnrichedGenesTable;
