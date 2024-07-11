import { DispatchContext } from "src/views/DifferentialExpression/common/store";
import { useContext, useEffect, useMemo, useState } from "react";

import { Props } from "./types";
import {
  setSelectedOptionsGroup1,
  setSelectedOptionsGroup2,
} from "src/views/DifferentialExpression/common/store/actions";

export const useConnect = ({
  options,
  allAvailableOptions,
  selectedOptionIds,
  queryGroupKey,
  isQueryGroup1,
}: {
  options: Props["options"];
  selectedOptionIds: Props["selectedOptionIds"];
  allAvailableOptions: Props["allAvailableOptions"];
  queryGroupKey: Props["queryGroupKey"];
  isQueryGroup1: Props["isQueryGroup1"];
}) => {
  const dispatch = useContext(DispatchContext);

  const [previousSelectedOptions, setPreviousSelectedOptions] = useState<
    Props["options"]
  >([]);

  const selectedOptions = useMemo(() => {
    const optionsMap = new Map(options.map((option) => [option.id, option]));
    const newOptions = allAvailableOptions
      .filter((option) => selectedOptionIds.includes(option.id))
      .map((availableOption) => {
        const foundOption = optionsMap.get(availableOption.id);
        return (
          foundOption || {
            id: availableOption.id,
            name: availableOption.name,
            unavailable: true,
          }
        );
      });
    newOptions.sort(
      (a, b) =>
        selectedOptionIds.indexOf(a.id) - selectedOptionIds.indexOf(b.id)
    );

    return newOptions;
  }, [options, allAvailableOptions, selectedOptionIds]);

  useEffect(() => {
    if (dispatch) {
      if (isQueryGroup1) {
        dispatch(setSelectedOptionsGroup1(queryGroupKey, selectedOptions));
      } else {
        dispatch(setSelectedOptionsGroup2(queryGroupKey, selectedOptions));
      }
    }
  }, [selectedOptions, dispatch, isQueryGroup1, queryGroupKey]);

  return {
    selectedOptions,
    previousSelectedOptions,
    setPreviousSelectedOptions,
  };
};
