import { useContext, useMemo } from "react";
import { StateContext } from "src/views/WheresMyGeneV2/common/store";
import {
  aggregateCollectionsFromDatasets,
  useFilterDimensions,
} from "src/common/queries/wheresMyGene";
import { Collections } from "./types";

export const useConnect = () => {
  const { selectedFilters } = useContext(StateContext);
  const { data: filterDimensions } = useFilterDimensions();

  const { datasets = [] } = filterDimensions;

  const collections: Collections = useMemo(() => {
    return aggregateCollectionsFromDatasets(datasets, selectedFilters);
  }, [datasets, selectedFilters]);

  return { collections };
};
