import React, { ReactElement, useState, useMemo, useEffect } from "react";
import {
  TableTitle,
  TableTitleWrapper,
  TableTitleInnerWrapper,
  TableUnavailableContainer,
  TableUnavailableHeader,
  TableUnavailableDescription,
} from "../common/style";
import { PublicationLinkWrapper } from "./style";
import Table from "../common/Table";
import Link from "../common/Link";
import { useCanonicalMarkers } from "src/common/queries/cellCards";
import DropdownSelect from "../common/DropdownSelect";
import { SelectChangeEvent } from "@mui/material";
import { Tooltip } from "@czi-sds/components";
import HelpTooltip from "../common/HelpTooltip";
import { ROUTES } from "src/common/constants/routes";

export const CELL_CARD_CANONICAL_MARKER_GENES_TABLE =
  "cell-card-canonical-marker-genes-table";
export const CELL_CARD_CANONICAL_MARKER_GENES_TABLE_DROPDOWN =
  "cell-card-canonical-marker-genes-table-dropdown";

interface TableRow {
  symbol: string;
  name: string;
  references: ReactElement | string;
}
const tableColumns: Array<keyof TableRow> = ["symbol", "name", "references"];

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
    // All Tissues selecteed by default
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

  const handleChange = (event: SelectChangeEvent<unknown>) => {
    setSelectedOrgan(event.target.value as string);
  };

  useEffect(() => {
    return () => {
      setSelectedOrgan("All Tissues");
    };
  }, []);
  const genesForShareUrl = tableRows.map((row) => row.symbol).join("%2C");

  return (
    <div data-testid={CELL_CARD_CANONICAL_MARKER_GENES_TABLE}>
      <TableTitleWrapper>
        <TableTitleInnerWrapper columnGap={4}>
          <TableTitle>Marker Genes</TableTitle>
          <HelpTooltip
            text={
              <div>
                {
                  "The below marker genes and associated publications were derived from the "
                }
                <Link
                  label={
                    "Anatomical Structures, Cell Types and Biomarkers (ASCT+B)"
                  }
                  url={"https://humanatlas.io/asctb-tables"}
                />
                {
                  " tables. The tables are authored and reviewed by an international team of anatomists, pathologists, physicians, and other experts."
                }
                <br />
                <br />
                <i>
                  Quardokus, Ellen, Bruce W. Herr II, Lisel Record, Katy Börner.
                  2022. HuBMAP ASCT+B Tables. Accessed May 16, 2023.
                </i>
              </div>
            }
          />
        </TableTitleInnerWrapper>
        {tableRows.length > 0 && (
          <TableTitleInnerWrapper>
            <DropdownSelect
              handleChange={handleChange}
              options={uniqueOrgans}
              selectedOption={selectedOrgan}
              testId={CELL_CARD_CANONICAL_MARKER_GENES_TABLE_DROPDOWN}
            />

            <Link
              url={`${ROUTES.WHERE_IS_MY_GENE}?genes=${genesForShareUrl}&ver=2`}
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
            Canonical marker genes for this cell type are unavailable at this
            time.
          </TableUnavailableDescription>
        </TableUnavailableContainer>
      )}
    </div>
  );
};
export default CanonicalMarkerGeneTable;
