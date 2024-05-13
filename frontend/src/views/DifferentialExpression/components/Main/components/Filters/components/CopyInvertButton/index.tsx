import { useContext, useEffect, useState } from "react";
import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";
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
import { StyledButtonIcon } from "./style";

interface Props {
  queryGroupKey: keyof QueryGroup;
  isCopy: boolean;
}
function CopyInvertButton({ queryGroupKey, isCopy }: Props): JSX.Element {
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
    if (isCopy) {
      const options = queryGroup1[queryGroupKey];
      const optionsWithNames = availableOptionsGroup2.filter((option) =>
        options.includes(option.id)
      );
      setOptions(optionsWithNames);
    } else {
      const options = queryGroup1[queryGroupKey];
      const optionsWithNames = availableOptionsGroup2.filter(
        (option) => !options.includes(option.id)
      );
      setOptions(optionsWithNames);
    }
  }, [availableOptionsGroup2, isCopy, queryGroup1, queryGroupKey]);

  useEffect(() => {
    setDisabled(options.length === 0);
  }, [options]);
  const handleClick = () => {
    dispatch && dispatch(selectQueryGroup2Filters(queryGroupKey, options));
  };

  return (
    <StyledButtonIcon
      onClick={handleClick}
      sdsIcon={isCopy ? "copy" : "circlesOverlap"}
      sdsSize="small"
      disabled={disabled}
    />
  );
}
export default CopyInvertButton;
