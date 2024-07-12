import { useMemo } from "react";
import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";
import { Tooltip } from "@czi-sds/components";
import { StyledTag } from "./style";
import { QUERY_GROUP_KEY_TO_TAG_SUFFIX_MAP } from "./constants";
import { getCellTypeLink } from "src/views/CellGuide/common/utils";
import { NO_ORGAN_ID } from "src/views/CellGuide/components/CellGuideCard/components/MarkerGeneTables/constants";
import Link from "src/views/CellGuide/components/CellGuideCard/components/common/Link";

const QueryGroupTags = ({
  submittedQueryGroup,
  submittedQueryGroupWithNames,
}: {
  submittedQueryGroup: QueryGroup;
  submittedQueryGroupWithNames: QueryGroup;
}) => {
  const nonEmptyQueryGroupKeys = useMemo(() => {
    return Object.keys(submittedQueryGroup).filter(
      (key) => submittedQueryGroup[key as keyof QueryGroup].length > 0
    );
  }, [submittedQueryGroup]);

  return (
    <>
      {nonEmptyQueryGroupKeys.map((key) => {
        const queryGroupKey = key as keyof QueryGroup;
        const selectedIds = submittedQueryGroup[queryGroupKey];
        const selectedNames = submittedQueryGroupWithNames[queryGroupKey];

        const suffix = QUERY_GROUP_KEY_TO_TAG_SUFFIX_MAP[queryGroupKey];
        const label =
          selectedNames.length > 1
            ? `${selectedNames.length} ${suffix}`
            : selectedNames[0];

        const getValue = (index: number) => {
          return key === "cellTypes" ? (
            <Link
              label={selectedNames[index]}
              url={getCellTypeLink({
                cellTypeId: selectedIds[index],
                tissueId: NO_ORGAN_ID,
              })}
            />
          ) : (
            selectedNames[index]
          );
        };
        const clickToViewText = "Click to view in CellGuide";
        const tooltipContent =
          selectedNames.length === 1 && key === "cellTypes" ? (
            clickToViewText
          ) : (
            <div>
              {key === "cellTypes" && <b>{clickToViewText}</b>}
              {selectedNames.map((value, index) => (
                <div key={`value-${value}-${index}`}>{getValue(index)}</div>
              ))}
            </div>
          );

        const isSingleCellType =
          key === "cellTypes" && selectedNames.length === 1;
        return (
          <Tooltip
            key={`${key}-tooltip`}
            sdsStyle="light"
            placement="top"
            width="wide"
            leaveDelay={0}
            disableHoverListener={
              key !== "cellTypes" && selectedNames.length === 1
            }
            title={tooltipContent}
          >
            <span>
              <TagWrapper
                key={`${key}-tag-wrapper`}
                selectedId={selectedIds[0]}
                isSingleCellType={isSingleCellType}
              >
                <StyledTag
                  color="neutral"
                  sdsType="secondary"
                  isSingleCellType={isSingleCellType}
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
