import { useContext } from "react";
import { QueryGroup } from "../../common/store/reducer";
import { StateContext } from "../../common/store";

export const useCraftPayloadWithQueryGroups = (): Record<string, string> => {
  const {
    selectedOptionsGroup1,
    selectedOptionsGroup2,
    excludeOverlappingCells,
  } = useContext(StateContext);

  const payload: Record<string, string> = {};

  Object.keys(selectedOptionsGroup1).forEach((key: string) => {
    const options1 = selectedOptionsGroup1[key as keyof QueryGroup]
      .filter((item) => !item.unavailable)
      .map((item) => item.id)
      .sort();
    const options2 = selectedOptionsGroup2[key as keyof QueryGroup]
      .filter((item) => !item.unavailable)
      .map((item) => item.id)
      .sort();
    if (options1.length > 0 || options2.length > 0) {
      payload[key] = `CG1: ${options1.join(",")}; CG2: ${options2.join(",")}`;
    }
  });

  payload["overlap"] = excludeOverlappingCells;
  return payload;
};
