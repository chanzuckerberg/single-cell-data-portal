import { useState, useEffect, useContext, useMemo, useCallback } from "react";
import { useDifferentialExpression } from "src/common/queries/differentialExpression";
import { StateContext } from "src/views/DifferentialExpression/common/store";
import {
  StyledIcon,
  ButtonsWrapper,
  ButtonLabel,
  InstructionsBody,
  InstructionsHeader,
  InstructionsWrapper,
  ResultsHeader,
  ResultsWrapper,
  ResultsHeaderWrapper,
  FlexRow,
} from "./style";

import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";
import { StyledSidebarDrawer } from "src/views/WheresMyGeneV2/components/Main/style";
import { DrawerSize } from "@blueprintjs/core";
import SourceData from "./components/SourceData";
import { DifferentialExpressionRow } from "./types";
import DifferentialExpressionResults from "./components/DifferentialExpressionResults";
import { applyFilter } from "./utils";
import {
  DIFFERENTIAL_EXPRESSION_INSTRUCTIONS_SIDEBAR,
  DIFFERENTIAL_EXPRESSION_RESULTS_DOWNLOAD_BUTTON,
  DIFFERENTIAL_EXPRESSION_SOURCE_DATA_BUTTON,
} from "src/views/DifferentialExpression/common/constants";

interface DeResultsProps {
  setIsLoading: (isLoading: boolean) => void;
}
export default function DeResults({
  setIsLoading,
}: DeResultsProps): JSX.Element {
  const [searchQuery, setSearchQuery] = useState("");
  const [lfcFilter, setLfcFilter] = useState("");
  const [effectSizeFilter, setEffectSizeFilter] = useState("");
  const [sortDirection, setSortDirection] = useState<"asc" | "desc">("desc");

  const [isSourceDatasetSidebarOpen, setIsSourceDatasetSidebarOpen] =
    useState(false);
  const { data, isLoading } = useDifferentialExpression();
  const {
    differentialExpressionResults: rawDifferentialExpressionResults,
    nOverlap,
  } = data;

  const [differentialExpressionResults, setDifferentialExpressionResults] =
    useState<DifferentialExpressionRow[]>([]);
  const {
    organismId,
    submittedQueryGroups: queryGroups,
    submittedQueryGroupsWithNames: queryGroupsWithNames,
  } = useContext(StateContext);

  useEffect(() => {
    if (!organismId || isLoading) return;

    // map ids to name
    const formattedDeResults = rawDifferentialExpressionResults.map(
      (diffExpResult) => {
        return {
          name: diffExpResult.gene_symbol,
          logFoldChange: diffExpResult.log_fold_change.toFixed(3),
          effectSize: diffExpResult.effect_size.toFixed(3),
          id: diffExpResult.gene_ontology_term_id,
          adjustedPValue: diffExpResult.adjusted_p_value.toExponential(3),
        };
      }
    );

    setDifferentialExpressionResults(formattedDeResults);
  }, [rawDifferentialExpressionResults, isLoading, organismId]);

  const namesToShow: string[][] = [];
  const { queryGroup1, queryGroup2 } = queryGroupsWithNames ?? {};
  for (const [index, queryGroupWithNames] of [
    queryGroup1,
    queryGroup2,
  ].entries()) {
    namesToShow.push([]);
    for (const key in queryGroupWithNames) {
      for (const value of queryGroupWithNames[key as keyof QueryGroup]) {
        namesToShow[index].push(value);
      }
    }
  }

  const showEmpty = !queryGroups;

  useEffect(() => {
    setIsLoading(isLoading);
  }, [isLoading, setIsLoading]);

  const sortedAndFilteredResults = useMemo(() => {
    return differentialExpressionResults
      .filter((result) =>
        searchQuery
          .toLowerCase()
          .split(",")
          .some((query) => result.name.toLowerCase().includes(query.trim()))
      )
      .filter((result) =>
        applyFilter(parseFloat(result.logFoldChange), lfcFilter)
      )
      .filter((result) =>
        applyFilter(parseFloat(result.effectSize), effectSizeFilter)
      )
      .sort((a, b) => {
        if (sortDirection === "asc") {
          return parseFloat(a.effectSize) - parseFloat(b.effectSize);
        } else if (sortDirection === "desc") {
          return parseFloat(b.effectSize) - parseFloat(a.effectSize);
        }
        return 0;
      });
  }, [
    differentialExpressionResults,
    searchQuery,
    lfcFilter,
    effectSizeFilter,
    sortDirection,
  ]);

  const downloadCSV = useCallback(() => {
    if (!queryGroupsWithNames) return;
    const { queryGroup1, queryGroup2 } = queryGroupsWithNames;

    const formatQueryGroupFilters = (queryGroup: QueryGroup): string => {
      return Object.entries(queryGroup)
        .filter(([_, values]) => values.length > 0)
        .map(([key, values]) => `${key}: ${values.join(", ")}`)
        .join(" | ");
    };

    const queryGroup1Filters = formatQueryGroupFilters(queryGroup1);
    const queryGroup2Filters = formatQueryGroupFilters(queryGroup2);

    const headers = [
      "Gene",
      "Log Fold Change",
      "Effect Size",
      "Adjusted P-Value",
    ];
    const csvRows = [
      `# Query Group 1 Filters: ${queryGroup1Filters}`,
      `# Query Group 2 Filters: ${queryGroup2Filters}`,
      headers.join(","),
    ];

    sortedAndFilteredResults.forEach((row) => {
      const values = [
        row.name,
        row.logFoldChange,
        row.effectSize,
        row.adjustedPValue,
      ];
      csvRows.push(values.join(","));
    });

    const csvContent = csvRows.join("\n");
    const blob = new Blob([csvContent], { type: "text/csv;charset=utf-8;" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.setAttribute("href", url);
    link.setAttribute("download", "differential_expression_results.csv");
    link.style.visibility = "hidden";
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  }, [sortedAndFilteredResults, queryGroupsWithNames]);

  return (
    <div>
      {!showEmpty ? (
        <ResultsWrapper>
          <ResultsHeaderWrapper>
            <ResultsHeader>Results</ResultsHeader>
            <FlexRow>
              <ButtonsWrapper
                onClick={downloadCSV}
                data-testid={DIFFERENTIAL_EXPRESSION_RESULTS_DOWNLOAD_BUTTON}
              >
                <StyledIcon sdsIcon="download" sdsSize="l" sdsType="static" />
                <ButtonLabel>Download</ButtonLabel>
              </ButtonsWrapper>
              <ButtonsWrapper
                onClick={() => setIsSourceDatasetSidebarOpen(true)}
                data-testid={DIFFERENTIAL_EXPRESSION_SOURCE_DATA_BUTTON}
              >
                <StyledIcon sdsIcon="infoCircle" sdsSize="l" sdsType="static" />
                <ButtonLabel>Source Data</ButtonLabel>
              </ButtonsWrapper>
            </FlexRow>
          </ResultsHeaderWrapper>
          {!!queryGroups &&
            !!queryGroupsWithNames &&
            !!organismId &&
            !isLoading && (
              <DifferentialExpressionResults
                queryGroups={queryGroups}
                queryGroupsWithNames={queryGroupsWithNames}
                organismId={organismId}
                sortedAndFilteredResults={sortedAndFilteredResults}
                nCellsOverlap={nOverlap}
                setSearchQuery={setSearchQuery}
                setLfcFilter={setLfcFilter}
                setEffectSizeFilter={setEffectSizeFilter}
                sortDirection={sortDirection}
                setSortDirection={setSortDirection}
              />
            )}
        </ResultsWrapper>
      ) : (
        <InstructionsWrapper
          data-testid={DIFFERENTIAL_EXPRESSION_INSTRUCTIONS_SIDEBAR}
        >
          <InstructionsHeader>Instructions</InstructionsHeader>
          <InstructionsBody>
            <ol>
              <li>
                Select a cell group of interest within the Cell Group 1 box by
                using the dropdown selectors.
                <br />
                <br />
                To copy the same selection over to Cell Group 2, click the copy
                button to the right of each dropdown in Cell Group 1.
              </li>
              <li>
                Within Cell Group 2, select a group that the cell group of
                interest will be compared to.
              </li>
            </ol>
          </InstructionsBody>
        </InstructionsWrapper>
      )}

      <StyledSidebarDrawer
        position="right"
        isOpen={isSourceDatasetSidebarOpen}
        title="Source Data"
        canEscapeKeyClose={true}
        canOutsideClickClose={true}
        onClose={() => setIsSourceDatasetSidebarOpen(false)}
        size={DrawerSize.SMALL}
      >
        <SourceData />
      </StyledSidebarDrawer>
    </div>
  );
}
