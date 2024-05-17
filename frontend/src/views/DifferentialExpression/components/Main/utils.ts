import { QueryGroups, QueryGroup } from "../../common/store/reducer";

export const craftPayloadWithQueryGroups = (
  queryGroups: QueryGroups
): Record<string, string> => {
  const payload: Record<string, string> = {};

  Object.keys(queryGroups.queryGroup1).forEach((key: string) => {
    payload[key] = `CG1: ${queryGroups.queryGroup1[
      key as keyof QueryGroup
    ].join(",")}, CG2: ${queryGroups.queryGroup2[key as keyof QueryGroup].join(
      ","
    )}`;
  });

  return payload;
};
