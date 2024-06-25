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
              {key === "cellTypes" && <b>{clickToViewText}</b>}
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
              <TagWrapper
                key={`${key}-tag-wrapper`}
                selectedId={selectedId[0]}
                isSingleCellType={key === "cellTypes" && selected.length === 1}
              >
                <StyledTag
                  color="neutral"
                  sdsType="secondary"
                  isSingleCellType={
                    key === "cellTypes" && selected.length === 1
                  }
                  label={label}
                />
              </TagWrapper>
            </span>
          </Tooltip>
        );
      })}
    </>
  );
};

interface TagWrapperProps {
  children: React.ReactNode;
  key: string;
  selectedId: string;
  isSingleCellType: boolean;
}

function TagWrapper({
  children,
  key,
  selectedId,
  isSingleCellType,
}: TagWrapperProps) {
  if (!isSingleCellType) {
    return <span key={key}>{children}</span>;
  }

  const url = getCellTypeLink({
    cellTypeId: selectedId,
    tissueId: NO_ORGAN_ID,
  });

  return (
    <a key={key} href={url} rel="noopener" target="_blank">
      {children}
    </a>
  );
}

export default QueryGroupTags;
