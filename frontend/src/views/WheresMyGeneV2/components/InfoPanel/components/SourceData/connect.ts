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
  let { datasets = [] } = filterDimensions;
  if (selectedFilters.datasets.length > 0)
    datasets = datasets.filter((dataset) =>
      selectedFilters.datasets.includes(dataset.id)
    );
  const collections: Collections = useMemo(() => {
    return aggregateCollectionsFromDatasets(datasets);
  }, [datasets]);
  return { collections };
};
