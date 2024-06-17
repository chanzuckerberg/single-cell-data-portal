import { useState, useMemo } from "react";

import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";

import { generateAndCopyShareUrl } from "./utils";

import { MAX_NUM_TOP_GENES_TO_PORT_TO_GE, ROWS_PER_PAGE } from "./constants";

import useProcessedQueryGroupFilterDimensions from "../../../common/query_group_filter_dimensions";
import { Props } from "./types";
import { useQueryGroupFilterDimensions } from "src/common/queries/differentialExpression";
import { EMPTY_ARRAY } from "src/common/constants/utils";

export const useConnect = ({
  queryGroups,
  queryGroupsWithNames,
  organismId,
  sortedAndFilteredResults,
  nCellsOverlap,
  sortDirection,
}: {
  queryGroups: Props["queryGroups"];
  queryGroupsWithNames: Props["queryGroupsWithNames"];
  organismId: Props["organismId"];
  sortedAndFilteredResults: Props["sortedAndFilteredResults"];
  nCellsOverlap: Props["nCellsOverlap"];
  sortDirection: Props["sortDirection"];
}) => {
  const [page, setPage] = useState(1);

  const { n_cells: nCellsGroup1 } = useProcessedQueryGroupFilterDimensions(
    queryGroups.queryGroup1
  );
  const { n_cells: nCellsGroup2 } = useProcessedQueryGroupFilterDimensions(
    queryGroups.queryGroup2
  );
  const handlePageChange = (
    _event: React.ChangeEvent<unknown>,
    page: number
  ) => {
    setPage(page);
  };

  const [openInGeHref1, openInGeHref2] = useMemo(() => {
    const generateUrl = (
      queryGroup: QueryGroup,
      queryGroupWithNames: QueryGroup,
      sliceFromBeginning: boolean
    ) => {
      return organismId
        ? generateAndCopyShareUrl({
            queryGroup: queryGroup,
            organism: organismId,
            genes: sliceFromBeginning
              ? sortedAndFilteredResults
                  .slice(0, MAX_NUM_TOP_GENES_TO_PORT_TO_GE)
                  .map((row) => row.name)
              : sortedAndFilteredResults
                  .slice(-MAX_NUM_TOP_GENES_TO_PORT_TO_GE)
                  .reverse()
                  .map((row) => row.name),
            cellTypes: queryGroupWithNames.cellTypes,
          })
        : "";
    };
    return [
      generateUrl(
        queryGroups.queryGroup1,
        queryGroupsWithNames.queryGroup1,
        sortDirection === "desc"
      ),
      generateUrl(
        queryGroups.queryGroup2,
        queryGroupsWithNames.queryGroup2,
        sortDirection === "asc"
      ),
    ];
  }, [
    queryGroups,
    queryGroupsWithNames,
    organismId,
    sortDirection,
    sortedAndFilteredResults,
  ]);

  const pageCount = Math.ceil(sortedAndFilteredResults.length / ROWS_PER_PAGE);

  const overlapPercent = (
    (nCellsOverlap / Math.max(nCellsGroup1, nCellsGroup2)) *
    100
  ).toFixed(2);

  const { data: filterDimensions1 } = useQueryGroupFilterDimensions(
    queryGroups.queryGroup1
  );
  const { data: filterDimensions2 } = useQueryGroupFilterDimensions(
    queryGroups.queryGroup2
  );

  const { datasets: datasets1 = EMPTY_ARRAY } = filterDimensions1;
  const { datasets: datasets2 = EMPTY_ARRAY } = filterDimensions2;
  return {
    page,
    setPage,
    sortDirection,
    nCellsGroup1,
    nCellsGroup2,
    openInGeHref1,
    openInGeHref2,
    pageCount,
    handlePageChange,
    overlapPercent,
    numDatasets1: `${datasets1.length} dataset${
      datasets1.length !== 1 ? "s" : ""
    }`,
    numDatasets2: `${datasets2.length} dataset${
      datasets2.length !== 1 ? "s" : ""
    }`,
  };
};
