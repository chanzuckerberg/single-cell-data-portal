import { useMemo, useState } from "react";

import { Props } from "./types";

export const useConnect = ({
  options,
  allAvailableOptions,
  selectedOptionIds,
}: {
  options: Props["options"];
  selectedOptionIds: Props["selectedOptionIds"];
  allAvailableOptions: Props["allAvailableOptions"];
}) => {
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
  return {
    selectedOptions,
    previousSelectedOptions,
    setPreviousSelectedOptions,
  };
};
