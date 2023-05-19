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
import { Tooltip } from "czifui";
import HelpTooltip from "../common/HelpTooltip";

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
    }
    if (!selectedOrgan) setSelectedOrgan("All Tissues");
    const uniqueOrgans = Array.from(organs);
    // put "All Tissues" first in the array
    return uniqueOrgans.sort((a, b) => {
      if (a === "All Tissues") return -1;
      if (b === "All Tissues") return 1;
      return a.localeCompare(b);
    });
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
      setSelectedOrgan("All Tissues");
    };
  }, []);
  // const genesForShareUrl = tableRows.map((row) => row.symbol).join("%2C");

  return (
    <div>
      <TableTitleWrapper>
        <TableTitleInnerWrapper columnGap={4}>
          <TableTitle id={MARKER_GENES_SECTION_ID}>Marker Genes</TableTitle>
          <HelpTooltip
            text={
              <>
                {
                  "The below marker genes and associated publications were derived from the "
                }
                <Link
                  label={
                    "Anatomical Structures, Cell Types and Biomarkers (ASCT+B)"
                  }
                  url={
                    "https://hubmapconsortium.github.io/ccf/pages/ccf-anatomical-structures.html"
                  }
                />
                {
                  " tables. The tables are authored and reviewed by an international team of anatomists, pathologists, physicians, and other experts."
                }
              </>
            }
          />
        </TableTitleInnerWrapper>
        {tableRows.length > 0 && (
          <TableTitleInnerWrapper>
            <DropdownSelect
              handleChange={handleChange}
              options={uniqueOrgans}
              selectedOption={selectedOrgan}
            />

            <Link
              url={`https://cellxgene.cziscience.com/gene-expression`}
              label={"Open in Gene Expression"}
            />
          </TableTitleInnerWrapper>
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
