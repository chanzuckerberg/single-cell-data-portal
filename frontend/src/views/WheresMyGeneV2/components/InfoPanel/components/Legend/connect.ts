import { useContext } from "react";
import { StateContext } from "src/views/WheresMyGeneV2/common/store";
import { useFilterDimensions } from "src/common/queries/wheresMyGene";

export const useConnect = () => {
  const { selectedFilters } = useContext(StateContext);
  const { data: filterDimensions } = useFilterDimensions();

  let { datasets = [] } = filterDimensions;

  if (selectedFilters.datasets.length > 0) {
    datasets = datasets.filter((dataset) =>
      selectedFilters.datasets.includes(dataset.id)
    );
  }

  const referenceCount = datasets.length;

  return { referenceCount };
};
