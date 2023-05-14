import React, { ReactElement, useMemo } from "react";
import {
  TableTitle,
  TableTitleWrapper,
  PublicationLinkWrapper,
  WmgLink,
  TableUnavailableContainer,
  TableUnavailableHeader,
  TableUnavailableDescription,
} from "../common/style";
import Table from "../Table";
import { MARKER_GENES_SECTION_ID } from "../CellCardSidebar";
import Link from "../common/Link";
import { useCanonicalMarkers } from "src/common/queries/cellCards";

interface TableRow {
  symbol: ReactElement;
  name: string;
  publications: ReactElement | string;
  organ: string;
  tissue: string;
}
const tableColumns: Array<keyof TableRow> = [
  "symbol",
  "name",
  "organ",
  "tissue",
  "publications",
];

interface Props {
  cellTypeId: string;
}

const CanonicalMarkerGeneTable = ({ cellTypeId }: Props) => {
  const { data: genes } = useCanonicalMarkers(cellTypeId);
  const tableRows: TableRow[] = useMemo(() => {
    if (!genes) return [];
    const rows = [];

    for (const markerGene of genes) {
      const publications = markerGene.publication.split(";;");
      const publicationTitles = markerGene.publication_titles.split(";;");
      const publicationLinks = (
        <PublicationLinkWrapper>
          {publications.map((publication, index) => {
            if (publication && publicationTitles[index]) {
              return (
                <Link
                  key={`${publication}-${index}`}
                  title={publicationTitles[index]}
                  url={`https://doi.org/${publication}`}
                />
              );
            }
          })}
        </PublicationLinkWrapper>
      );

      rows.push({
        symbol: (
          <Link
            title={`${markerGene.symbol}`}
            url={`https://www.genecards.org/cgi-bin/carddisp.pl?gene=${markerGene.symbol}`}
          />
        ),
        name: markerGene.name,
        publications: publicationLinks,
        organ: markerGene.tissue_general,
        tissue: markerGene.tissue_specific,
      });
    }
    return rows;
  }, [genes]);

  const genesForShareUrl = tableRows.map((row) => row.symbol).join("%2C");

  return (
    <div>
      <TableTitleWrapper>
        <TableTitle id={MARKER_GENES_SECTION_ID}>Marker Genes</TableTitle>
        {tableRows.length > 0 && (
          <WmgLink
            href={`https://cellxgene.cziscience.com/gene-expression?tissues=lung&genes=${genesForShareUrl}&ver=2`}
            target="_blank"
          >
            Open in Gene Expression
          </WmgLink>
        )}
      </TableTitleWrapper>
      {tableRows.length ? (
        <Table<TableRow> columns={tableColumns} rows={tableRows} />
      ) : (
        <TableUnavailableContainer>
          <TableUnavailableHeader>
            No canonical marker genes
          </TableUnavailableHeader>
          <TableUnavailableDescription>
            Canonical marker genes genes for this cell type are unavailable at
            this time.
          </TableUnavailableDescription>
        </TableUnavailableContainer>
      )}
    </div>
  );
};
export default CanonicalMarkerGeneTable;
