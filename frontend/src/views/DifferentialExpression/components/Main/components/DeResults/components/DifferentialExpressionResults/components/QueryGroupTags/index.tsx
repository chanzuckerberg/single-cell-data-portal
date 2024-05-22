import { useMemo } from "react";
import {
  QueryGroup,
  QueryGroups,
} from "src/views/DifferentialExpression/common/store/reducer";
import { Tooltip } from "@czi-sds/components";
import { StyledTag } from "./style";
import { QUERY_GROUP_KEY_TO_TAG_SUFFIX_MAP } from "./constants";

const QueryGroupTags = ({
  queryGroupsWithNames,
  isQueryGroup1,
}: {
  queryGroupsWithNames: QueryGroups;
  isQueryGroup1?: boolean;
}) => {
  const queryGroup = isQueryGroup1
    ? queryGroupsWithNames.queryGroup1
    : queryGroupsWithNames.queryGroup2;
  const nonEmptyQueryGroupKeys = useMemo(() => {
    return Object.keys(queryGroup).filter(
      (key) => queryGroup[key as keyof QueryGroup].length > 0
    );
  }, [queryGroup]);

  return (
    <>
      {nonEmptyQueryGroupKeys.map((key) => {
        const queryGroupKey = key as keyof QueryGroup;
        const selected = queryGroup[queryGroupKey];
        const suffix = QUERY_GROUP_KEY_TO_TAG_SUFFIX_MAP[queryGroupKey];
        const label = `${selected.length} ${
          selected.length > 1 ? suffix.plural : suffix.single
        }`;
        return (
          <Tooltip
            key={`${key}-tooltip`}
            sdsStyle="light"
            placement="top"
            width="wide"
            leaveDelay={0}
            title={selected.map((value, index) => (
              <div key={`value-${value}-${index}`}>{value}</div>
            ))}
          >
            <span>
              <StyledTag
                key={key}
                color="gray"
                sdsType="secondary"
                label={label}
              />
            </span>
          </Tooltip>
        );
      })}
    </>
  );
};

export default QueryGroupTags;
