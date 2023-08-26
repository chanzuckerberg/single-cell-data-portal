import { useMemo } from "react";
import {
  SourceCollectionsQueryResponse,
  SourceCollectionsQueryResponseEntry,
} from "src/common/queries/cellGuide";
import { ALL_TISSUES } from "../../MarkerGeneTables/constants";
import { filterDescendantsOfAncestorTissueId } from "src/views/CellGuide/common/utils";

function _passSelectionCriteria({
  collection,
  selectedOrganismLabel,
  selectedOrganLabel,
  selectedOrganId,
}: {
  collection: SourceCollectionsQueryResponseEntry;
  selectedOrganismLabel: string;
  selectedOrganLabel: string;
  selectedOrganId: string;
}): boolean {
  const collectionOrganismLabels = collection.organism.map(
    (organism) => organism.label
  );
  const collectionTissueIds = collection.tissue.map(
    (tissue) => tissue.ontology_term_id
  );

  if (!collectionOrganismLabels.includes(selectedOrganismLabel)) return false;

  if (selectedOrganLabel === ALL_TISSUES) {
    return true;
  }

  const descendantTissueIds = filterDescendantsOfAncestorTissueId(
    collectionTissueIds,
    selectedOrganId
  );

  return descendantTissueIds.length > 0;
}

function _filterCollections({
  collections,
  selectedOrganismLabel,
  selectedOrganLabel,
  selectedOrganId,
}: {
  collections: SourceCollectionsQueryResponse;
  selectedOrganismLabel: string;
  selectedOrganLabel: string;
  selectedOrganId: string;
}): SourceCollectionsQueryResponse {
  const filteredCollections: SourceCollectionsQueryResponse = [];

  for (const collection of collections) {
    if (
      _passSelectionCriteria({
        collection: collection,
        selectedOrganismLabel: selectedOrganismLabel,
        selectedOrganLabel: selectedOrganLabel,
        selectedOrganId: selectedOrganId,
      })
    )
      filteredCollections.push(collection);
  }

  return filteredCollections;
}

function _sortCollections(
  collections: SourceCollectionsQueryResponse
): SourceCollectionsQueryResponse {
  return collections.sort((a, b) => {
    const aOrganisms = a.organism.map((organism) => organism.label);
    const bOrganisms = b.organism.map((organism) => organism.label);
    if (aOrganisms.length === 0 && bOrganisms.length === 0) return 0;
    if (aOrganisms.includes("Homo sapiens")) return -1;
    if (bOrganisms.includes("Homo sapiens")) return 1;
    const aFirstOrganism = aOrganisms.at(0) ?? "";
    const bFirstOrganism = bOrganisms.at(0) ?? "";
    return aFirstOrganism.localeCompare(bFirstOrganism);
  });
}

export function useDataSourceFilter({
  collections,
  selectedOrganismLabel,
  selectedOrganId,
  selectedOrganLabel,
}: {
  collections: SourceCollectionsQueryResponse;
  selectedOrganismLabel: string;
  selectedOrganId: string;
  selectedOrganLabel: string;
}): SourceCollectionsQueryResponseEntry[] {
  return useMemo(() => {
    if (!collections) return [];

    const filteredCollections = _filterCollections({
      collections,
      selectedOrganismLabel,
      selectedOrganLabel,
      selectedOrganId,
    });
    return _sortCollections(filteredCollections);
  }, [collections, selectedOrganismLabel, selectedOrganLabel, selectedOrganId]);
}
