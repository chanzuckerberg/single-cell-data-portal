import React, { ReactElement } from "react";
import { WmgLink } from "./style";
import {
  TableTitle,
  TableTitleWrapper,
  PublicationLinkWrapper,
} from "../common/style";
import { allCellTypeMarkerGenes } from "src/views/CellCards/common/fixtures";
import Table from "../Table";

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
  id: string;
}

const CanonicalMarkerGeneTable = ({ cellTypeId, id }: Props) => {
  const tableRows: TableRow[] = [];
  if (cellTypeId in allCellTypeMarkerGenes) {
    const genes =
      allCellTypeMarkerGenes[cellTypeId as keyof typeof allCellTypeMarkerGenes];
    for (const markerGene of genes) {
      const publications = markerGene.publication.split(";;");
      const publicationTitles = markerGene.publication_titles.split(";;");
      const publicationLinks = (
        <PublicationLinkWrapper>
          {publications.map((publication, index) => {
            if (publication && publicationTitles[index]) {
              return (
                <Link
                  key={publication}
                  title={publicationTitles[index]}
                  url={`https://doi.org/${publication}`}
                />
              );
            }
          })}
        </PublicationLinkWrapper>
      );

      tableRows.push({
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
  }

  const genesForShareUrl = tableRows.map((row) => row.symbol).join("%2C");

  return (
    <div id={id}>
      <TableTitleWrapper>
        <TableTitle>Marker Genes</TableTitle>
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
        <div>Canonical marker genes are not available yet.</div>
      )}
    </div>
  );
};
export default CanonicalMarkerGeneTable;
