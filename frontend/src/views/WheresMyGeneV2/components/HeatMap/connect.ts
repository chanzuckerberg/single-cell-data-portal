import {
  useCallback,
  useContext,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";
import { EMPTY_ARRAY, EMPTY_OBJECT } from "src/common/constants/utils";
import {
  CellTypeRow,
  OntologyTerm,
  generateTermsByKey,
  usePrimaryFilterDimensions,
} from "src/common/queries/wheresMyGene";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGeneV2/common/store";
import {
  useSortedGeneNames,
  useTissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap,
} from "src/views/WheresMyGeneV2/components/HeatMap/hooks/useSortedGeneNames";
import {
  useHandleExpandedTissueIds,
  useTrackHeatMapLoaded,
} from "src/views/WheresMyGeneV2/components/HeatMap/hooks";

import { Props } from "./types";
import { DefaultAutocompleteOption } from "@czi-sds/components";
import {
  CellType,
  GeneExpressionSummary,
  SORT_BY,
  Tissue,
} from "src/views/WheresMyGeneV2/common/types";
import {
  addCellInfoCellType,
  setFilteredCellTypes,
  toggleExpandedTissueId,
} from "src/views/WheresMyGeneV2/common/store/actions";
import { useSortedCellTypesByTissueName } from "src/views/WheresMyGeneV2/components/HeatMap/hooks/useSortedCellTypesByTissueName";
import { cloneDeep } from "lodash";

export function useConnect({
  cellTypes,
  cellTypeSortBy,
  genes,
  geneSortBy,
  selectedGeneExpressionSummariesByTissueName,
  setTissuesByName,
  tissuesByName,
}: Props) {
  const {
    expandedTissueIds,
    filteredCellTypeIds,
    filteredCellTypes,
    selectedFilters: { tissues: filteredTissueIds },
    xAxisHeight,
  } = useContext(StateContext);

  const selectedCellTypeOptions = useMemo(() => {
    return filteredCellTypes.map((cellType) => ({
      name: cellType,
    }));
  }, [filteredCellTypes]);

  // Loading state per tissue
  const [isLoading, setIsLoading] = useState(setInitialIsLoading(cellTypes));
  const chartWrapperRef = useRef<HTMLDivElement>(null);
  const dispatch = useContext(DispatchContext);

  const { data } = usePrimaryFilterDimensions(); //temp explicit version

  // Get tissueName to ID map for use in find marker genes
  useEffect(() => {
    let result: { [name: string]: OntologyTerm } = EMPTY_OBJECT;

    if (data) {
      const { tissues } = data;

      result = generateTermsByKey(tissues, "name");
    }

    setTissuesByName(result);
  }, [data, setTissuesByName]);

  const cellTypesByName = useMemo(() => {
    const result: { [name: string]: CellType } = {};

    Object.values(cellTypes).forEach((cellTypes) => {
      cellTypes.forEach((cellType) => {
        result[cellType.cellTypeName] = cellType;
      });
    });

    return result;
  }, [cellTypes]);

  const generateMarkerGenes = (cellType: CellType, tissueID: string) => {
    dispatch?.(addCellInfoCellType({ cellType, tissueID }));
  };

  const tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap =
    useTissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap(
      selectedGeneExpressionSummariesByTissueName
    );

  const sortedGeneNames = useSortedGeneNames({
    genes,
    geneSortBy,
    selectedCellTypes: cellTypes,
    tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap,
  });

  const sortedCellTypesByTissueName = useSortedCellTypesByTissueName({
    cellTypeSortBy,
    genes,
    selectedCellTypes: cellTypes,
    tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap,
  });

  const geneNameToIndex = useMemo(() => {
    const result: { [key: string]: number } = {};

    for (const [index, gene] of Object.entries(sortedGeneNames)) {
      result[gene] = Number(index);
    }

    return result;
  }, [sortedGeneNames]);

  const orderedSelectedGeneExpressionSummariesByTissueName = useMemo(() => {
    const result: { [tissueName: string]: GeneExpressionSummary[] } = {};

    for (const [tissueName, geneExpressionSummary] of Object.entries(
      selectedGeneExpressionSummariesByTissueName
    )) {
      // (thuang): sort() mutates the array, so we need to clone it
      result[tissueName] = cloneDeep(
        geneExpressionSummary.sort((a, b) => {
          if (!a || !b) return -1;

          return geneNameToIndex[a.name] - geneNameToIndex[b.name];
        })
      );
    }

    return result;
  }, [selectedGeneExpressionSummariesByTissueName, geneNameToIndex]);

  /**
   * (thuang): Tissues to display after applying filters
   */
  const displayedTissues = useMemo(() => {
    return Object.values(tissuesByName)
      .filter(({ id }) => {
        return !filteredTissueIds.length || filteredTissueIds.includes(id);
      })
      .filter(({ name }) => {
        if (!filteredCellTypes.length) return true;

        const tissueCellTypes = sortedCellTypesByTissueName[name];

        return tissueCellTypes?.some((cellType) => {
          return filteredCellTypes.includes(cellType.cellTypeName);
        });
      });
  }, [
    filteredTissueIds,
    filteredCellTypes,
    sortedCellTypesByTissueName,
    tissuesByName,
  ]);

  const displayedTissueIds = useMemo(() => {
    return displayedTissues.map(({ id }) => id);
  }, [displayedTissues]);

  /**
   * (thuang): Derive displayed cell types from `displayedTissues`,
   * `expandedTissueIds`, and `filteredCellTypes`
   */
  const displayedCellTypes = useMemo(() => {
    const result = new Set<string>();

    displayedTissues.forEach(({ id, name }) => {
      result.add(id + id);

      if (expandedTissueIds.includes(id)) {
        const tissueCellTypes = sortedCellTypesByTissueName[name];

        tissueCellTypes?.forEach((cellType) => {
          if (
            !filteredCellTypes.length ||
            filteredCellTypes.includes(cellType.cellTypeName)
          ) {
            result.add(id + cellType.cellTypeName);
          }
        });
      }
    });

    return result;
  }, [
    displayedTissues,
    expandedTissueIds,
    filteredCellTypes,
    sortedCellTypesByTissueName,
  ]);

  const handleExpandCollapse = useCallback(
    (tissueId: string, tissueName: Tissue) => {
      dispatch?.(toggleExpandedTissueId({ tissueId, tissueName }));
    },
    [dispatch]
  );

  const uniqueCellTypes = useMemo(() => {
    const result: Set<string> = new Set<string>();
    Object.values(sortedCellTypesByTissueName).forEach((cellTypes) => {
      cellTypes.forEach((cellType) => {
        if (!cellType.cellTypeName.includes("UBERON:"))
          result.add(cellType.cellTypeName);
      });
    });
    return [...result].sort().map((cellType) => ({ name: cellType }));
  }, [sortedCellTypesByTissueName]);

  const handleFilteredCellTypesChange = useCallback(
    (_: unknown, rawNewFilteredCellTypes: DefaultAutocompleteOption[]) => {
      if (!Object.keys(cellTypesByName).length) return;

      const newCellTypeNames = rawNewFilteredCellTypes.map(
        (cellType) => cellType.name
      );
      const cellTypeIds = newCellTypeNames.map(
        (name) => cellTypesByName[name].id
      );

      /**
       * (thuang): Don't dispatch if the new filtered cell types are the same.
       * Otherwise, it will cause infinite loop.
       */
      if (String(filteredCellTypes) === String(newCellTypeNames)) return;

      dispatch?.(
        setFilteredCellTypes({
          filteredCellTypes: newCellTypeNames,
          filteredCellTypeIds: cellTypeIds,
          displayedTissueIds,
        })
      );
    },
    [cellTypesByName, displayedTissueIds, dispatch, filteredCellTypes]
  );

  const handleCellTypeDelete = (cellTypeNameToDelete: string) => () => {
    const cellTypeIdToDelete = cellTypesByName[cellTypeNameToDelete].id;
    const newCellTypeNames = filteredCellTypes.filter(
      (cellType) => !(cellTypeNameToDelete === cellType)
    );
    const newCellTypeIds = filteredCellTypeIds.filter(
      (cellTypeId) => !(cellTypeIdToDelete === cellTypeId)
    );

    dispatch?.(
      setFilteredCellTypes({
        filteredCellTypes: newCellTypeNames,
        filteredCellTypeIds: newCellTypeIds,
        displayedTissueIds,
      })
    );
  };

  /**
   * All tissue cell types to render in YAxisCharts
   */
  const allTissueCellTypes = useMemo(() => {
    return displayedTissues
      .sort((a, b) => {
        // sort tissues alphabetically
        return a.name.localeCompare(b.name);
      })
      .flatMap((tissue: OntologyTerm) => {
        const { id, name } = tissue;

        const tissueCellTypes = getTissueCellTypes({
          cellTypeSortBy,
          cellTypes,
          sortedCellTypesByTissueName,
          tissue: name,
          tissueID: id,
          displayedCellTypes,
        });

        return tissueCellTypes.length > 0
          ? {
              tissueId: id,
              tissueName: name,
              tissueCellTypes,
            }
          : [];
      });
  }, [
    cellTypeSortBy,
    cellTypes,
    sortedCellTypesByTissueName,
    displayedCellTypes,
    displayedTissues,
  ]);

  const totalElementsCount = useMemo(() => {
    return allTissueCellTypes.reduce((count, element) => {
      if (element.tissueCellTypes.length > 0) {
        return count + element.tissueCellTypes.length + 1;
      }
      return count++;
    }, 0);
  }, [allTissueCellTypes]);

  useTrackHeatMapLoaded({
    selectedGenes: genes,
    displayedCellTypes,
    selectedCellTypes: filteredCellTypes,
  });

  useHandleExpandedTissueIds({
    filteredCellTypeIds,
    filteredTissueIds,
    displayedTissueIds,
    dispatch,
  });

  return {
    allTissueCellTypes,
    chartWrapperRef,
    expandedTissueIds,
    filteredCellTypes,
    generateMarkerGenes,
    handleCellTypeDelete,
    handleExpandCollapse,
    handleFilteredCellTypesChange,
    isLoading,
    orderedSelectedGeneExpressionSummariesByTissueName,
    selectedCellTypeOptions,
    setIsLoading,
    sortedGeneNames,
    totalElementsCount,
    uniqueCellTypes,
    useHandleExpandedTissueIds,
    useTrackHeatMapLoaded,
    xAxisHeight,
  };
}

function setInitialIsLoading(cellTypes: Props["cellTypes"]) {
  return Object.keys(cellTypes).reduce((isLoading, tissue) => {
    return { ...isLoading, [tissue]: false };
  }, {});
}

function getTissueCellTypes({
  cellTypes,
  sortedCellTypesByTissueName,
  tissue,
  tissueID,
  cellTypeSortBy,
  displayedCellTypes,
}: {
  cellTypes: { [tissue: Tissue]: CellTypeRow[] };
  sortedCellTypesByTissueName: { [tissue: string]: CellTypeRow[] };
  tissue: Tissue;
  tissueID: string;
  cellTypeSortBy: SORT_BY;
  displayedCellTypes: Set<string>;
}) {
  const tissueCellTypes = cellTypes[tissue];

  if (!tissueCellTypes || tissueCellTypes.length === 0) return EMPTY_ARRAY;

  const sortedTissueCellTypes = sortedCellTypesByTissueName[tissue];

  let ret =
    (cellTypeSortBy === SORT_BY.CELL_ONTOLOGY
      ? tissueCellTypes
      : sortedTissueCellTypes) || EMPTY_ARRAY;

  ret = ret.filter(
    (cellType) =>
      displayedCellTypes.has(tissueID + cellType.cellTypeName) &&
      cellType.total_count > 9
  );

  return ret;
}
