import { useMemo } from "react";
import { CanonicalMarkersQueryResponse } from "src/common/queries/cellGuide";
import { HOMO_SAPIENS } from "../constants";

interface CanonicalMarkerGeneTableData {
  symbol: string;
  name: string;
  referenceData: {
    publicationTitles: string[];
    publicationTitlesToIndex: Map<string, number>;
    publications: string[];
    numReferences: number;
  };
}

function _getSortedOrgans(genes: CanonicalMarkersQueryResponse): string[] {
  const organs = new Set<string>();
  for (const markerGene of genes) {
    organs.add(markerGene.tissue);
  }
  return Array.from(organs).sort((a, b) => {
    return a.localeCompare(b);
  });
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

export function useCanonicalMarkerGenesTableRowsAndFilters({
  genes,
  selectedOrgan,
}: {
  genes: CanonicalMarkersQueryResponse;
  selectedOrgan: string;
}): {
  canonicalMarkerGeneTableData: CanonicalMarkerGeneTableData[];
  uniqueOrganisms: string[];
} {
  return useMemo(() => {
    if (!genes)
      return {
        canonicalMarkerGeneTableData: [],
        uniqueOrganisms: [],
      };

    // get sorted organs
    const sortedOrgans = _getSortedOrgans(genes);
    const selectedOrganFilter = !sortedOrgans.includes(selectedOrgan)
      ? ""
      : selectedOrgan;

    const sortedOrganisms = [HOMO_SAPIENS]; // ASCTB tables only report data for humans

    const rows: CanonicalMarkerGeneTableData[] = [];
    const publicationTitlesToIndex = _getPublicationTitlesToIndex(
      genes,
      selectedOrganFilter
    );

    for (const markerGene of genes) {
      const { tissue, publication, publication_titles, symbol, name } =
        markerGene;

      if (tissue !== selectedOrganFilter) continue;
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
      uniqueOrganisms: sortedOrganisms,
    };
  }, [genes, selectedOrgan]);
}
