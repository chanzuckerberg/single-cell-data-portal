import { useContext, useEffect, useState } from "react";
import {
  DispatchContext,
  StateContext,
} from "src/views/DifferentialExpression/common/store";

import {
  submitQueryGroups,
  clearQueryGroup1Filters,
  clearQueryGroup2Filters,
} from "src/views/DifferentialExpression/common/store/actions";
import useProcessedQueryGroupFilterDimensions from "./components/common/query_group_filter_dimensions";

import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { craftPayloadWithQueryGroups } from "./utils";

export const useConnect = () => {
  const [isLoading, setIsLoading] = useState<boolean>(false);
  const [isLoadingGetDeQuery, setIsLoadingGetDeQuery] =
    useState<boolean>(false);

  useEffect(() => {
    setIsLoadingGetDeQuery(isLoadingGetDeQuery);
  }, [isLoadingGetDeQuery]);
  const dispatch = useContext(DispatchContext);
  const { queryGroups } = useContext(StateContext);
  const { queryGroup1, queryGroup2 } = queryGroups;

  // check if any values in queryGroup1 are not empty
  const isQueryGroup1NotEmpty = Object.values(queryGroup1).some(
    (value) => value.length > 0
  );
  const isQueryGroup2NotEmpty = Object.values(queryGroup2).some(
    (value) => value.length > 0
  );
  const canRunDifferentialExpression =
    !isLoading && isQueryGroup1NotEmpty && isQueryGroup2NotEmpty;

  const handleRunDifferentialExpression = () => {
    if (!dispatch) return;
    dispatch(submitQueryGroups());

    track(
      EVENTS.DE_FIND_GENES_CLICKED,
      craftPayloadWithQueryGroups(queryGroups)
    );
  };

  const handleClearQueryGroups = () => {
    if (!dispatch) return;
    dispatch(clearQueryGroup1Filters());
    dispatch(clearQueryGroup2Filters());
  };

  const { n_cells: nCellsGroup1, isLoading: isLoadingGroup1 } =
    useProcessedQueryGroupFilterDimensions(queryGroup1);
  const { n_cells: nCellsGroup2, isLoading: isLoadingGroup2 } =
    useProcessedQueryGroupFilterDimensions(queryGroup2);

  return {
    isLoading,
    setIsLoading,
    queryGroup1,
    queryGroup2,
    canRunDifferentialExpression,
    handleRunDifferentialExpression,
    handleClearQueryGroups,
    nCellsGroup1,
    isLoadingGroup1,
    nCellsGroup2,
    isLoadingGroup2,
  };
};
