import { useContext, useEffect, useState } from "react";
import useProcessedQueryGroupFilterDimensions, {
  FilterOptionDimensions,
} from "../../../common/query_group_filter_dimensions";
import {
  DispatchContext,
  StateContext,
} from "src/views/DifferentialExpression/common/store";
import { selectQueryGroup2Filters } from "src/views/DifferentialExpression/common/store/actions";
import { QUERY_GROUP_KEY_TO_FILTER_DIMENSION_MAP } from "../../../common/constants";
import { FilterOption } from "../../types";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { Props } from "./types";

export const useConnect = ({
  queryGroupKey,
}: {
  queryGroupKey: Props["queryGroupKey"];
}) => {
  const {
    queryGroups: { queryGroup1, queryGroup2 },
  } = useContext(StateContext);
  const dispatch = useContext(DispatchContext);
  const [disabled, setDisabled] = useState(false);
  const [options, setOptions] = useState<FilterOption[]>([]);
  const filterDimensionKey = QUERY_GROUP_KEY_TO_FILTER_DIMENSION_MAP[
    queryGroupKey
  ] as keyof FilterOptionDimensions;

  const {
    availableFilters: { [filterDimensionKey]: availableOptionsGroup2 },
  } = useProcessedQueryGroupFilterDimensions(queryGroup2);

  useEffect(() => {
    const options = queryGroup1[queryGroupKey];
    const optionsWithNames = availableOptionsGroup2.filter((option) =>
      options.includes(option.id)
    );
    setOptions(optionsWithNames);
  }, [availableOptionsGroup2, queryGroup1, queryGroupKey]);

  useEffect(() => {
    setDisabled(options.length === 0);
  }, [options]);
  const handleClick = () => {
    dispatch && dispatch(selectQueryGroup2Filters(queryGroupKey, options));
    track(EVENTS.DE_CG_COPY_CLICKED, {
      category: queryGroupKey,
      values: options.map((o) => o.id).join(","),
    });
  };

  return {
    disabled,
    handleClick,
  };
};
