import React, { ReactElement, useMemo } from "react";
import { TableTitle, TableTitleWrapper, WmgLink } from "../common/style";
import Table from "../Table";
import Link from "../common/Link";
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

  const tableRows: TableRow[] = useMemo(() => {
    if (!genes) return [];
    const rows = [];
    for (const markerGene of genes) {
      const { pc, me, name, symbol } = markerGene;
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
  }, [genes]);

  const genesForShareUrl = tableRows.map((row) => row.symbol).join("%2C");

  return (
    <>
      <TableTitleWrapper>
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
