import { useCallback, useContext } from "react";

import { DispatchContext } from "src/views/DifferentialExpression/common/store";

import {
  EMPTY_FILTERS,
  QueryGroup,
} from "src/views/DifferentialExpression/common/store/reducer";
import {
  selectQueryGroup1Filters,
  selectQueryGroup2Filters,
} from "src/views/DifferentialExpression/common/store/actions";

import { FilterOption, Props } from "./types";
import useProcessedQueryGroupFilterDimensions from "../common/query_group_filter_dimensions";

import { QUERY_GROUP_KEYS_TO_FILTER_EVENT_MAP } from "./constants";
import { track } from "src/common/analytics";

export const useConnect = ({ queryGroup, isQueryGroup1 }: Props) => {
  const dispatch = useContext(DispatchContext);
  const selectQueryGroupFilters = isQueryGroup1
    ? selectQueryGroup1Filters
    : selectQueryGroup2Filters;

  const { availableFilters } =
    useProcessedQueryGroupFilterDimensions(queryGroup);

  const handleFilterChange = useCallback(
    function handleFilterChange_(
      key: Partial<keyof QueryGroup>
    ): (options: FilterOption[] | undefined) => void {
      return (options: FilterOption[] | undefined): void => {
        if (!dispatch || !options) {
          return;
        }
        const event = QUERY_GROUP_KEYS_TO_FILTER_EVENT_MAP[key];
        event &&
          options.length > 0 &&
          track(event, {
            [key]: options.map((o) => o.id).join(","),
          });
        dispatch(selectQueryGroupFilters(key, options));
      };
    },
    [dispatch, selectQueryGroupFilters]
  );

  const { availableFilters: allAvailableFilters } =
    useProcessedQueryGroupFilterDimensions(EMPTY_FILTERS);

  return {
    availableFilters,
    allAvailableFilters,
    handleFilterChange,
  };
};
