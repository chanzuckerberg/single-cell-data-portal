import React, { ReactElement, useState, useMemo, useEffect } from "react";
import {
  TableTitle,
  TableTitleWrapper,
  PublicationLinkWrapper,
  TableUnavailableContainer,
  TableUnavailableHeader,
  TableUnavailableDescription,
} from "../common/style";
import Table from "../Table";
import { MARKER_GENES_SECTION_ID } from "../CellCardSidebar";
import Link from "../common/Link";
import { useCanonicalMarkers } from "src/common/queries/cellCards";
import { TableTitleInnerWrapper } from "../EnrichedGenesTable/style";
import DropdownSelect from "../common/DropdownSelect";
import { SelectChangeEvent } from "@mui/material";

interface TableRow {
  symbol: string;
  name: string;
  publications: ReactElement | string;
}
const tableColumns: Array<keyof TableRow> = ["symbol", "name", "publications"];

interface Props {
  cellTypeId: string;
}

const CanonicalMarkerGeneTable = ({ cellTypeId }: Props) => {
  const { data: genes } = useCanonicalMarkers(cellTypeId);
  const [selectedOrgan, setSelectedOrgan] = useState("");

  const uniqueOrgans = useMemo(() => {
    if (!genes) return [];
    const organs = new Set<string>();
    for (const markerGene of genes) {
      organs.add(markerGene.tissue_general);
      if (!selectedOrgan) setSelectedOrgan(markerGene.tissue_general);
    }
    return Array.from(organs);
  }, [genes, cellTypeId]);

  const tableRows: TableRow[] = useMemo(() => {
    if (!genes) return [];
    const rows = [];

    for (const markerGene of genes) {
      if (markerGene.tissue_general !== selectedOrgan) continue;
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
        symbol: markerGene.symbol,
        name: markerGene.name,
        publications: publicationLinks,
      });
    }
    return rows;
  }, [genes, selectedOrgan, cellTypeId]);

  const handleChange = (event: SelectChangeEvent) => {
    setSelectedOrgan(event.target.value as string);
  };

  useEffect(() => {
    return () => {
      setSelectedOrgan("");
    };
  }, []);
  // const genesForShareUrl = tableRows.map((row) => row.symbol).join("%2C");

  return (
    <div>
      <TableTitleWrapper>
        <TableTitleInnerWrapper>
          <TableTitle id={MARKER_GENES_SECTION_ID}>Marker Genes</TableTitle>
          {tableRows.length > 0 && (
            <DropdownSelect
              handleChange={handleChange}
              options={uniqueOrgans}
              selectedOption={selectedOrgan}
            />
          )}
        </TableTitleInnerWrapper>
        {tableRows.length > 0 && (
          <Link
            url={`https://cellxgene.cziscience.com/gene-expression`}
            title={"Open in Gene Expression"}
          />
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
