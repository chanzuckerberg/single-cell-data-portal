import { useMemo } from "react";
import {
  CellGroupTitle,
  CellGroupTitleWrapper,
  CellGroupWrapper,
  EffectSizeHeaderWrapper,
  EffectSizeIndicator,
  FilterTagsWrapper,
  StyledTextField,
  TableHeaderWrapper,
  TableWrapper,
  OpenInGE,
  StyledCallout,
  StyledIconImage,
  StyledTooltipText,
  StyledInterpretButton,
  GroupButtonsWrapper,
  Overlay,
} from "./style";
import cxgIcon from "./images/cxg.svg";
import { Pagination } from "@mui/material";
import Table from "src/views/CellGuide/components/CellGuideCard/components/common/Table";
import { ButtonIcon, Tooltip } from "@czi-sds/components";

import { DifferentialExpressionRow } from "../../types";
import { MAX_NUM_TOP_GENES_TO_PORT_TO_GE, ROWS_PER_PAGE } from "./constants";
import { CellCountTitle, Spinner } from "../../../../style";
import QueryGroupTags from "./components/QueryGroupTags";
import {
  DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_INFO,
  DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_INFO,
  DIFFERENTIAL_EXPRESSION_EFFECT_SIZE_FILTER,
  DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT,
  DIFFERENTIAL_EXPRESSION_GENES_FILTER,
  DIFFERENTIAL_EXPRESSION_LFC_FILTER,
  DIFFERENTIAL_EXPRESSION_OPEN_IN_GE_1_BUTTON,
  DIFFERENTIAL_EXPRESSION_OPEN_IN_GE_2_BUTTON,
  DIFFERENTIAL_EXPRESSION_RESULTS_CALLOUT,
  DIFFERENTIAL_EXPRESSION_RESULTS_TABLE,
  DIFFERENTIAL_EXPRESSION_SORT_DIRECTION,
} from "src/views/DifferentialExpression/common/constants";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { useConnect } from "./connect";
import { Props } from "./types";
import InterpretationCard from "./components/InterpretationCard";

const DifferentialExpressionResults = ({
  queryGroups,
  queryGroupsWithNames,
  organismId,
  sortedAndFilteredResults,
  nCellsOverlap,
  setSearchQuery,
  setLfcFilter,
  setEffectSizeFilter,
  sortDirection,
  setSortDirection,
}: Props) => {
  const {
    page,
    setPage,
    nCellsGroup1,
    nCellsGroup2,
    openInGeHref1,
    openInGeHref2,
    pageCount,
    handlePageChange,
    overlapPercent,
    genesToInterpret,
    isQueryGroup1BeingInterpreted,
    setIsQueryGroup1BeingInterpreted,
    isLoadingInterpret,
    setIsLoadingInterpret,
    interpretationCardVisible,
    setInterpretationCardVisible,
  } = useConnect({
    queryGroups,
    queryGroupsWithNames,
    organismId,
    sortedAndFilteredResults,
    nCellsOverlap,
    sortDirection,
  });

  const columnIdToName: Record<
    keyof Omit<DifferentialExpressionRow, "adjustedPValue">,
    string | JSX.Element
  > = useMemo(() => {
    const handleSearch = (event: React.ChangeEvent<HTMLInputElement>) => {
      setSearchQuery(event.target.value);
      track(EVENTS.DE_SEARCH_GENE, { gene: event.target.value });
      setPage(1);
    };
    const handleLfcFilter = (event: React.ChangeEvent<HTMLInputElement>) => {
      setLfcFilter(event.target.value);
      setPage(1);
    };
    const handleEffectSizeFilter = (
      event: React.ChangeEvent<HTMLInputElement>
    ) => {
      setEffectSizeFilter(event.target.value);
      setPage(1);
    };

    const handleSortDirectionChange = () => {
      setSortDirection((prevDirection) =>
        prevDirection === "asc" ? "desc" : "asc"
      );
    };

    return {
      name: (
        <TableHeaderWrapper data-testid={DIFFERENTIAL_EXPRESSION_GENES_FILTER}>
          Gene
          <StyledTextField
            variant="outlined"
            onChange={handleSearch}
            placeholder="e.g. JCHAIN"
          />
        </TableHeaderWrapper>
      ),
      logFoldChange: (
        <TableHeaderWrapper data-testid={DIFFERENTIAL_EXPRESSION_LFC_FILTER}>
          Log Fold Change
          <StyledTextField
            placeholder="e.g >1.0"
            variant="outlined"
            onChange={handleLfcFilter}
          />
        </TableHeaderWrapper>
      ),
      effectSize: (
        <TableHeaderWrapper
          data-testid={DIFFERENTIAL_EXPRESSION_EFFECT_SIZE_FILTER}
        >
          <EffectSizeHeaderWrapper
            data-testid={DIFFERENTIAL_EXPRESSION_SORT_DIRECTION}
            onClick={handleSortDirectionChange}
          >
            Effect Size{" "}
            <ButtonIcon
              sdsIcon={sortDirection === "asc" ? "chevronUp" : "chevronDown"}
              sdsSize="small"
            />
          </EffectSizeHeaderWrapper>
          <StyledTextField
            placeholder="e.g >1.0"
            variant="outlined"
            onChange={handleEffectSizeFilter}
          />
        </TableHeaderWrapper>
      ),
    };
  }, [
    sortDirection,
    setSearchQuery,
    setLfcFilter,
    setEffectSizeFilter,
    setSortDirection,
    setPage,
  ]);

  const interpretText1 =
    isQueryGroup1BeingInterpreted && isLoadingInterpret
      ? "Interpreting..."
      : "Interpret";

  const interpretText2 =
    !isQueryGroup1BeingInterpreted && isLoadingInterpret
      ? "Interpreting..."
      : "Interpret";

  return (
    <>
      {interpretationCardVisible && !isLoadingInterpret && <Overlay />}
      <CellGroupWrapper data-testid={DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_INFO}>
        <CellGroupTitleWrapper>
          <CellGroupTitle>Cell Group 1</CellGroupTitle>
          <GroupButtonsWrapper>
            <StyledInterpretButton
              onClick={() => {
                setIsQueryGroup1BeingInterpreted(true);
                setInterpretationCardVisible(true);
              }}
            >
              {isLoadingInterpret && isQueryGroup1BeingInterpreted && (
                <Spinner />
              )}
              {interpretText1}
            </StyledInterpretButton>
            <Tooltip
              sdsStyle="dark"
              placement="left"
              width="wide"
              arrow
              title={
                <StyledTooltipText>
                  Launch Gene Expression with the top{" "}
                  {MAX_NUM_TOP_GENES_TO_PORT_TO_GE} genes having the largest
                  positive effect sizes, indicating an increase in expression.
                  <br />
                  <br />
                  Filters applied to Cell Group 1 will also be used in Gene
                  Expression.
                </StyledTooltipText>
              }
            >
              <OpenInGE
                href={openInGeHref1}
                target="_blank"
                rel="noopener noreferrer"
                data-testid={DIFFERENTIAL_EXPRESSION_OPEN_IN_GE_1_BUTTON}
              >
                <StyledIconImage alt="CxG icon" src={cxgIcon} />
              </OpenInGE>
            </Tooltip>
          </GroupButtonsWrapper>
        </CellGroupTitleWrapper>
        <CellCountTitle data-testid={DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT}>
          {nCellsGroup1.toLocaleString()} cells |{" "}
          <EffectSizeIndicator>{"(+) Effect Size"}</EffectSizeIndicator>
        </CellCountTitle>
        <FilterTagsWrapper>
          <QueryGroupTags
            isQueryGroup1
            queryGroupsWithNames={queryGroupsWithNames}
          />
        </FilterTagsWrapper>
      </CellGroupWrapper>
      <CellGroupWrapper data-testid={DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_INFO}>
        <CellGroupTitleWrapper>
          <CellGroupTitle>Cell Group 2</CellGroupTitle>
          <GroupButtonsWrapper>
            <StyledInterpretButton
              sdsStyle="primary"
              sdsType="filled"
              onClick={() => {
                setIsQueryGroup1BeingInterpreted(false);
                setInterpretationCardVisible(true);
              }}
            >
              {isLoadingInterpret && !isQueryGroup1BeingInterpreted && (
                <Spinner />
              )}
              {interpretText2}
            </StyledInterpretButton>
            <Tooltip
              sdsStyle="dark"
              placement="left"
              width="wide"
              arrow
              title={
                <StyledTooltipText>
                  Launch Gene Expression with the top{" "}
                  {MAX_NUM_TOP_GENES_TO_PORT_TO_GE} genes having the largest
                  negative effect sizes, indicating a decrease in expression.
                  <br />
                  <br />
                  Filters applied to Cell Group 2 will also be used in Gene
                  Expression.
                </StyledTooltipText>
              }
            >
              <OpenInGE
                href={openInGeHref2}
                target="_blank"
                rel="noopener noreferrer"
                data-testid={DIFFERENTIAL_EXPRESSION_OPEN_IN_GE_2_BUTTON}
              >
                <StyledIconImage alt="CxG icon" src={cxgIcon} />
              </OpenInGE>
            </Tooltip>
          </GroupButtonsWrapper>
        </CellGroupTitleWrapper>
        <CellCountTitle data-testid={DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT}>
          {nCellsGroup2.toLocaleString()} cells |{" "}
          <EffectSizeIndicator>{"(-) Effect Size"}</EffectSizeIndicator>
        </CellCountTitle>
        <FilterTagsWrapper>
          <QueryGroupTags queryGroupsWithNames={queryGroupsWithNames} />
        </FilterTagsWrapper>
      </CellGroupWrapper>
      {nCellsOverlap > 0 && (
        <StyledCallout
          data-testid={DIFFERENTIAL_EXPRESSION_RESULTS_CALLOUT}
          intent={parseFloat(overlapPercent) > 25 ? "warning" : "info"}
        >
          {nCellsOverlap.toLocaleString()} overlapping cells ({overlapPercent}
          %) between groups. Selecting highly overlapping groups can result in
          underestimation of differences.
        </StyledCallout>
      )}
      <TableWrapper data-testid={DIFFERENTIAL_EXPRESSION_RESULTS_TABLE}>
        <Table<Omit<DifferentialExpressionRow, "adjustedPValue">>
          columns={["name", "logFoldChange", "effectSize"]}
          rows={sortedAndFilteredResults.slice(
            (page - 1) * ROWS_PER_PAGE,
            page * ROWS_PER_PAGE
          )}
          hoverable={false}
          columnIdToName={columnIdToName}
        />

        <Pagination count={pageCount} page={page} onChange={handlePageChange} />
      </TableWrapper>
      {interpretationCardVisible && (
        <InterpretationCard
          setIsVisible={setInterpretationCardVisible}
          isQueryGroup1={isQueryGroup1BeingInterpreted}
          differentialExpressionResults={genesToInterpret}
          setIsLoadingInterpret={setIsLoadingInterpret}
        />
      )}
    </>
  );
};

export default DifferentialExpressionResults;
