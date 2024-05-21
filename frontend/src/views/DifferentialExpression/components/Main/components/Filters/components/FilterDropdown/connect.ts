import { useMemo } from "react";

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
  const selectedOptions = useMemo(() => {
    const optionsMap = new Map(options.map((option) => [option.id, option]));
    return allAvailableOptions
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
  }, [options, allAvailableOptions, selectedOptionIds]);
  return {
    selectedOptions,
  };
};
