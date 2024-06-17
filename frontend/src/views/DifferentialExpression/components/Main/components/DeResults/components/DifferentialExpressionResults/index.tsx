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
} from "./style";
import cxgIcon from "./images/cxg.svg";
import { Pagination } from "@mui/material";
import Table from "src/views/CellGuide/components/CellGuideCard/components/common/Table";
import { Button, Tooltip } from "@czi-sds/components";

import { DifferentialExpressionRow } from "../../types";
import { MAX_NUM_TOP_GENES_TO_PORT_TO_GE, ROWS_PER_PAGE } from "./constants";
import { CellCountTitle } from "../../../../style";
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
import HelpTooltip from "src/views/CellGuide/components/CellGuideCard/components/common/HelpTooltip";
import { ROUTES } from "src/common/constants/routes";

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
            Effect Size
            <HelpTooltip
              title={"Effect Size"}
              dark
              text={
                <>
                  The effect size is the log fold change of the gene expression
                  between two cell groups normalized by the pooled standard
                  deviation, otherwise known as Cohen&apos;s d.
                  <br />
                  <br />
                  <div>
                    <a href={ROUTES.FMG_DOCS} rel="noopener" target="_blank">
                      Learn more about how Cohen&apos;s d is calculated.
                    </a>
                  </div>
                </>
              }
            />
            <Button
              sdsStyle="icon"
              icon={sortDirection === "asc" ? "ChevronUp" : "ChevronDown"}
              sdsSize="small"
              sdsType="tertiary"
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

  return (
    <>
      <CellGroupWrapper data-testid={DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_INFO}>
        <CellGroupTitleWrapper>
          <CellGroupTitle>Cell Group 1</CellGroupTitle>
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
          intent={parseFloat(overlapPercent) > 25 ? "notice" : "info"}
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
    </>
  );
};

export default DifferentialExpressionResults;
