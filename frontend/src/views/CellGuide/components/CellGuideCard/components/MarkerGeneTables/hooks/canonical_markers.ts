import { useMemo } from "react";
import {
  CanonicalMarkersQueryResponse,
  CanonicalMarkersQueryResponseEntry,
  useAllTissuesLookupTables,
  useCanonicalMarkers,
} from "src/common/queries/cellGuide";
import {
  ALL_TISSUES,
  HOMO_SAPIENS,
  NO_ORGAN_ID,
  TISSUE_AGNOSTIC,
} from "../constants";
import { isTissueIdDescendantOfAncestorTissueId } from "src/views/CellGuide/common/utils";

export interface CanonicalMarkerGeneTableData {
  symbol: string;
  name: string;
  referenceData: {
    publicationTitles: string[];
    publicationTitlesToIndex: Map<string, number>;
    publications: string[];
    numReferences: number;
  };
}

function _getPublicationTitlesToIndex(
  genes: CanonicalMarkersQueryResponse,
  selectedOrganFilter: string
): Map<string, number> {
  const publicationTitlesToIndex = new Map();
  let index = 0;
  for (const markerGene of genes) {
    if (markerGene.tissue !== selectedOrganFilter) continue;
    const publicationTitles = markerGene.publication_titles.split(";;");
    for (let i = 0; i < publicationTitles.length; i += 1) {
      if (
        !publicationTitlesToIndex.has(publicationTitles[i]) &&
        publicationTitles[i] !== ""
      ) {
        publicationTitlesToIndex.set(publicationTitles[i], index);
        index += 1;
      }
    }
  }
  return publicationTitlesToIndex;
}

function _getReferenceData(
  publication: string,
  publication_titles: string,
  publicationTitlesToIndex: Map<string, number>
): CanonicalMarkerGeneTableData["referenceData"] {
  let publications = Array.from(new Set(publication.split(";;")));
  let publicationTitles = Array.from(new Set(publication_titles.split(";;")));

  const sortedPublicationsAndTitles = publications
    .map((pub, i) => [pub, publicationTitles[i]])
    .sort(
      (a, b) =>
        (publicationTitlesToIndex.get(a[1]) ?? 0) -
        (publicationTitlesToIndex.get(b[1]) ?? 0)
    )
    .filter(
      (publicationTitle, index) => publicationTitle && publications[index]
    );

  publications = sortedPublicationsAndTitles.map((pub) => pub[0]);
  publicationTitles = sortedPublicationsAndTitles.map((pub) => pub[1]);
  return {
    publicationTitles,
    publicationTitlesToIndex,
    publications,
    numReferences: sortedPublicationsAndTitles.length,
  };
}

function _passSelectionCriteria({
  markerGene,
  allTissuesLabelToIdMap,
  selectedOrganId,
}: {
  markerGene: CanonicalMarkersQueryResponseEntry;
  allTissuesLabelToIdMap: Map<string, string>;
  selectedOrganId: string;
}): boolean {
  // There are marker genes tissues labeled as "All Tissues"
  // so select them only when selectedOrganLabel is "All Tissues"
  if (selectedOrganId === NO_ORGAN_ID) {
    return markerGene.tissue === ALL_TISSUES;
  }

  const tissue_id = allTissuesLabelToIdMap.get(markerGene.tissue);

  if (!tissue_id) {
    return false;
  }

  return isTissueIdDescendantOfAncestorTissueId(tissue_id, selectedOrganId);
}

export function useCanonicalMarkerGenesTableRowsAndFilters({
  cellTypeId,
  organName,
  organismName,
  organId,
}: {
  cellTypeId: string;
  organName: string;
  organismName: string;
  organId: string;
}): {
  canonicalMarkerGeneTableData: CanonicalMarkerGeneTableData[];
} {
  const { allTissuesLabelToIdLookup: allTissuesLabelToIdMap } =
    useAllTissuesLookupTables(cellTypeId);
  const { data: genes } = useCanonicalMarkers(cellTypeId);
  return useMemo(() => {
    if (!genes || organismName != HOMO_SAPIENS)
      return {
        canonicalMarkerGeneTableData: [],
      };

    const rows: CanonicalMarkerGeneTableData[] = [];

    const organLabel = TISSUE_AGNOSTIC ? ALL_TISSUES : organName;
    const publicationTitlesToIndex = _getPublicationTitlesToIndex(
      genes,
      organLabel
    );

    for (const markerGene of genes) {
      if (
        !_passSelectionCriteria({
          markerGene: markerGene,
          allTissuesLabelToIdMap: allTissuesLabelToIdMap,
          selectedOrganId: organId,
        })
      )
        continue;

      const { publication, publication_titles, symbol, name } = markerGene;
      const referenceData = _getReferenceData(
        publication,
        publication_titles,
        publicationTitlesToIndex
      );
      rows.push({
        name,
        symbol,
        referenceData,
      });
    }

    // Sort rows by number of references
    if (rows.length) {
      rows.sort((a, b) => {
        return b.referenceData.numReferences - a.referenceData.numReferences;
      });
    }

    return {
      canonicalMarkerGeneTableData: rows,
    };
  }, [genes, organismName, organName, organId, allTissuesLabelToIdMap]);
}
