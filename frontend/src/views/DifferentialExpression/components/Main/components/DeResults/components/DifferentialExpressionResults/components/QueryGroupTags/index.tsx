import { useMemo } from "react";
import {
  QueryGroup,
  QueryGroups,
} from "src/views/DifferentialExpression/common/store/reducer";
import { Tooltip } from "@czi-sds/components";
import { StyledTag } from "./style";
import { QUERY_GROUP_KEY_TO_TAG_SUFFIX_MAP } from "./constants";
import { getCellTypeLink } from "src/views/CellGuide/common/utils";
import { NO_ORGAN_ID } from "src/views/CellGuide/components/CellGuideCard/components/MarkerGeneTables/constants";
import Link from "src/views/CellGuide/components/CellGuideCard/components/common/Link";

const QueryGroupTags = ({
  queryGroups,
  queryGroupsWithNames,
  isQueryGroup1,
}: {
  queryGroups: QueryGroups;
  queryGroupsWithNames: QueryGroups;
  isQueryGroup1?: boolean;
}) => {
  const queryGroup = isQueryGroup1
    ? queryGroupsWithNames.queryGroup1
    : queryGroupsWithNames.queryGroup2;
  const queryGroupIds = isQueryGroup1
    ? queryGroups.queryGroup1
    : queryGroups.queryGroup2;
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
        const selectedId = queryGroupIds[queryGroupKey];
        const suffix = QUERY_GROUP_KEY_TO_TAG_SUFFIX_MAP[queryGroupKey];
        const label =
          selected.length > 1 ? `${selected.length} ${suffix}` : selected[0];

        const getValue = (index: number) => {
          return key === "cellTypes" ? (
            <Link
              label={selected[index]}
              url={getCellTypeLink({
                cellTypeId: selectedId[index],
                tissueId: NO_ORGAN_ID,
              })}
            />
          ) : (
            selected[index]
          );
        };
        const clickToViewText = "Click to view in CellGuide";
        const tooltipContent =
          selected.length === 1 && key === "cellTypes" ? (
            clickToViewText
          ) : (
            <div>
              <b>{clickToViewText}</b>
              {selected.map((value, index) => (
                <div key={`value-${value}-${index}`}>{getValue(index)}</div>
              ))}
            </div>
          );

        return (
          <Tooltip
            key={`${key}-tooltip`}
            sdsStyle="light"
            placement="top"
            width="wide"
            leaveDelay={0}
            disableHoverListener={key !== "cellTypes" && selected.length === 1}
            title={tooltipContent}
          >
            <span>
              <StyledTag
                key={key}
                color="neutral"
                sdsType="secondary"
                isSingleCellType={key === "cellTypes" && selected.length === 1}
                label={label}
                onClick={
                  key === "cellTypes" && selected.length === 1
                    ? () => {
                        const url = getCellTypeLink({
                          cellTypeId: selectedId[0],
                          tissueId: NO_ORGAN_ID,
                        });
                        window.open(url, "_blank");
                      }
                    : undefined
                }
              />
            </span>
          </Tooltip>
        );
      })}
    </>
  );
};

export default QueryGroupTags;
