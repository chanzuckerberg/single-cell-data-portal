import { QueryGroups, QueryGroup } from "../../common/store/reducer";

export const craftPayloadWithQueryGroups = (
  queryGroups: QueryGroups
): Record<string, string> => {
  const payload: Record<string, string> = {};

  Object.keys(queryGroups.queryGroup1).forEach((key: string) => {
    const options1 = queryGroups.queryGroup1[key as keyof QueryGroup];
    const options2 = queryGroups.queryGroup2[key as keyof QueryGroup];
    if (options1.length > 0 || options2.length > 0) {
      payload[key] = `CG1: ${options1.join(",")}, CG2: ${options2.join(",")}`;
    }
  });

  return payload;
};
