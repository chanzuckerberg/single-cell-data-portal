import { useState, useEffect, useContext, useMemo, useCallback } from "react";
import { useDifferentialExpression } from "src/common/queries/differentialExpression";
import { StateContext } from "src/views/DifferentialExpression/common/store";
import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";
import { DifferentialExpressionRow, Props } from "./types";
import { applyFilter } from "./utils";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { craftPayloadWithQueryGroups } from "../../utils";

export const useConnect = ({ setIsLoading }: Props) => {
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
    const formattedResults = rawDifferentialExpressionResults.map(
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
    setDifferentialExpressionResults(formattedResults);
  }, [rawDifferentialExpressionResults, isLoading, organismId]);

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
    if (!queryGroupsWithNames || !queryGroups || isLoading) return;

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

    track(EVENTS.DE_DOWNLOAD_CLICKED, craftPayloadWithQueryGroups(queryGroups));
  }, [sortedAndFilteredResults, queryGroups, queryGroupsWithNames, isLoading]);

  return {
    queryGroups,
    isLoading,
    queryGroupsWithNames,
    organismId,
    setSearchQuery,
    setLfcFilter,
    setEffectSizeFilter,
    sortDirection,
    setSortDirection,
    isSourceDatasetSidebarOpen,
    setIsSourceDatasetSidebarOpen,
    sortedAndFilteredResults,
    downloadCSV,
    showEmpty,
    nOverlap,
  };
};
