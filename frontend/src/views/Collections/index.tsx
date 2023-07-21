import Head from "next/head";
import React, { useEffect, useMemo } from "react";
import { Column, Filters, useFilters, useSortBy, useTable } from "react-table";
import { PLURALIZED_METADATA_LABEL } from "src/common/constants/metadata";
import { ROUTES } from "src/common/constants/routes";
import { useCategoryFilter } from "src/common/hooks/useCategoryFilter/useCategoryFilter";
import { useSessionStorage } from "src/common/hooks/useSessionStorage";
import { useFetchCollectionRows } from "src/common/queries/filter";
import { KEYS } from "src/common/sessionStorage/set";
import { useExplainTombstoned } from "src/components/Collections/common/utils";
import { CollectionsGrid } from "src/components/Collections/components/Grid/components/CollectionsGrid/style";
import {
  CATEGORY_FILTER_ID,
  CategoryView,
  CellPropsValue,
  CollectionRow,
  MultiPanelSelectedUIState,
  RowPropsValue,
} from "src/components/common/Filter/common/entities";
import { ontologyLabelCellAccessorFn } from "src/components/common/Filter/common/utils";
import {
  arraySortingFn,
  buildTableCountSummary,
} from "src/components/common/Grid/common/utils";
import DiseaseCell from "src/components/common/Grid/components/DiseaseCell";
import { GridHero } from "src/components/common/Grid/components/Hero";
import LinkCell from "src/components/common/Grid/components/LinkCell";
import NTagCell from "src/components/common/Grid/components/NTagCell";
import { Title } from "src/components/common/Grid/components/Title";
import CreateCollection from "src/components/CreateCollectionModal";
import SideBar from "src/components/common/SideBar";
import { CollectionsView as View } from "./style";
import { RightAlignCell } from "src/components/common/Grid/components/RightAlignCell";
import CountCell from "src/components/common/Grid/components/CountCell";
import {
  CATEGORY_FILTER_DENY_LIST,
  CATEGORY_FILTER_PARTITION_LIST,
  COLLECTION_CELL_COUNT,
  COLLECTION_CURATOR_NAME,
  COLLECTION_ID,
  COLLECTION_NAME,
  COLLECTION_RECENCY,
  COLLECTION_REVISED_BY,
  COLLECTION_STATUS,
  COLLECTION_SUMMARY_CITATION,
  COLLECTIONS_COLUMN_DENY_LIST,
  COLUMN_DENY_LIST,
  COLUMN_ID_RECENCY,
} from "src/views/Collections/common/constants";
import { ALIGNMENT } from "src/components/common/Grid/common/entities";
import StatusCell from "src/components/common/Grid/components/StatusCell";
import RevisionButton from "src/components/common/Grid/components/RevisionButton";
import CategoryFilters from "src/components/common/Filter/components/Filters";
import { useViewMode, VIEW_MODE } from "src/common/hooks/useViewMode";
import Loader from "src/components/common/Grid/components/Loader";

export default function Collections(): JSX.Element {
  const { mode, status } = useViewMode();

  // Pop toast if user has been redirected from a tombstoned collection.
  useExplainTombstoned();

  // Filterable collection datasets joined from datasets index and collections (or user-collections) index responses.
  const {
    isError,
    isSuccess,
    rows: collectionRows,
  } = useFetchCollectionRows(mode, status);

  // Show collections list when collections are successfully loaded.
  const shouldShowCollections = !isError && isSuccess;

  // Loading indicator for curator mode, when collections are not yet successfully loaded.
  const shouldShowLoader = mode === VIEW_MODE.CURATOR && !shouldShowCollections;

  // Column configuration backing table.
  const columnConfig: Column<CollectionRow>[] = useMemo(
    () => [
      {
        Cell: ({ row }: RowPropsValue<CollectionRow>) => {
          return (
            <Title>
              <LinkCell
                data-testid="collection-link"
                url={ROUTES.COLLECTION.replace(":id", row.values.id)}
                value={row.values.name}
              />
            </Title>
          );
        },
        Header: "Collections",
        accessor: COLLECTION_NAME,
        disableSortBy: true,
        showCountAndTotal: true,
      },
      // Viewable only in curator mode, required for filter.
      {
        Cell: ({ row }: RowPropsValue<CollectionRow>) => {
          return (
            <StatusCell
              revisionButton={<RevisionButton id={row.values.revisedBy} />}
              status={row.values.STATUS}
            />
          );
        },
        Header: "Status",
        accessor: COLLECTION_STATUS,
        disableSortBy: false,
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.STATUS,
        sortType: arraySortingFn,
      },
      // Viewable only in curator mode, required for filter.
      {
        Header: "Curator",
        accessor: COLLECTION_CURATOR_NAME,
        disableSortBy: false,
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.CURATOR_NAME,
        sortType: "alphanumeric",
      },
      // Viewable in collections mode, hidden in curator mode.
      {
        Cell: ({ row }: RowPropsValue<CollectionRow>) => {
          return <div>{row.values.summaryCitation || "No publication"}</div>;
        },
        Header: "Publication",
        accessor: COLLECTION_SUMMARY_CITATION,
        disableSortBy: true,
      },
      {
        Cell: ({ value }: CellPropsValue<string[]>) => (
          <NTagCell label={PLURALIZED_METADATA_LABEL.TISSUE} values={value} />
        ),
        Header: "Tissue",
        accessor: ontologyLabelCellAccessorFn("tissue"),
        disableSortBy: true,
        id: "tissue",
      },
      {
        Cell: ({ value }: CellPropsValue<string[]>) => (
          <DiseaseCell
            label={PLURALIZED_METADATA_LABEL.DISEASE}
            values={value}
          />
        ),
        Header: "Disease",
        accessor: ontologyLabelCellAccessorFn("disease"),
        disableSortBy: true,
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.DISEASE,
      },
      // Viewable only in curator mode, required for filter.
      {
        Cell: ({ value }: CellPropsValue<string[]>) => (
          <NTagCell label={PLURALIZED_METADATA_LABEL.ASSAY} values={value} />
        ),
        Header: "Assay",
        accessor: ontologyLabelCellAccessorFn("assay"),
        disableSortBy: true,
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.ASSAY,
      },
      {
        Cell: ({ value }: CellPropsValue<string[]>) => (
          <NTagCell label={PLURALIZED_METADATA_LABEL.ORGANISM} values={value} />
        ),
        Header: "Organism",
        accessor: ontologyLabelCellAccessorFn("organism"),
        disableSortBy: true,
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.ORGANISM,
      },
      // Viewable only in curator mode.
      {
        Cell: ({ value }: CellPropsValue<number | null>) => (
          <RightAlignCell>
            <CountCell cellCount={value || 0} />
          </RightAlignCell>
        ),
        Header: "Cells",
        accessor: COLLECTION_CELL_COUNT,
        alignment: ALIGNMENT.RIGHT,
        disableSortBy: true,
        filter: "between",
        id: CATEGORY_FILTER_ID.CELL_COUNT,
      },
      // Hidden, required for sorting
      {
        accessor: COLLECTION_RECENCY,
        id: COLUMN_ID_RECENCY,
      },
      // Hidden, required for accessing collection ID via row.values, for building link to collection detail page.
      {
        accessor: COLLECTION_ID,
      },
      // Hidden, required for filter.
      {
        accessor: "cellTypeCalculated",
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.CELL_TYPE_CALCULATED,
      },
      // Hidden, required for filter.
      {
        accessor: ontologyLabelCellAccessorFn("self_reported_ethnicity"),
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.SELF_REPORTED_ETHNICITY,
      },
      {
        accessor: "development_stage_ancestors",
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
      },
      // Hidden, required for filter.
      {
        accessor: "publicationAuthors",
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.PUBLICATION_AUTHORS,
      },
      // Hidden, required for filter.
      {
        accessor: "publicationDateValues",
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.PUBLICATION_DATE_VALUES,
      },
      // Hidden, required for accessing revised by via row.values, for building status.
      {
        accessor: COLLECTION_REVISED_BY,
      },
      // Hidden, required for filter.
      {
        accessor: ontologyLabelCellAccessorFn("sex"),
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.SEX,
      },
      // Hidden, required for filter.
      {
        accessor: "suspension_type",
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.SUSPENSION_TYPE,
      },
      // Hidden, required for filter.
      {
        accessor: "tissueCalculated",
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.TISSUE_CALCULATED,
      },
    ],
    []
  );

  // Handle initial filter state and save of filter state beyond component scope.
  const [initialFilters, storeFilters] = useSessionStorage<
    Filters<CollectionRow>
  >(KEYS.FILTER_COLLECTIONS, []);

  // Handle initial multi-panel filter UI state and save of multi-panel filter UI state beyond component scope.
  const [initialMultiPanelSelectedUIState, storeMultiPanelSelectedUIState] =
    useSessionStorage<MultiPanelSelectedUIState>(
      KEYS.FILTER_COLLECTIONS_SELECTED_UI,
      {}
    );

  // Table init
  const tableInstance = useTable<CollectionRow>(
    {
      columns: columnConfig,
      data: collectionRows,
      disableSortBy: false,
      initialState: {
        filters: initialFilters,
        hiddenColumns: COLLECTIONS_COLUMN_DENY_LIST,
        sortBy: [
          {
            desc: true,
            id: COLUMN_ID_RECENCY,
          },
        ],
      },
    },
    useFilters,
    useSortBy
  );

  // Determine the set of categories to display for the collections view.
  const categories = useMemo<Set<CATEGORY_FILTER_ID>>(() => {
    return Object.values(CATEGORY_FILTER_ID)
      .filter(
        (categoryFilterId: CATEGORY_FILTER_ID) =>
          !CATEGORY_FILTER_DENY_LIST[mode].includes(categoryFilterId)
      )
      .reduce((accum, categoryFilterId: CATEGORY_FILTER_ID) => {
        accum.add(categoryFilterId);
        return accum;
      }, new Set<CATEGORY_FILTER_ID>());
  }, [mode]);

  // Determine the hidden columns.
  const hiddenColumns = useMemo<string[]>(() => COLUMN_DENY_LIST[mode], [mode]);

  // Filter init.
  const {
    preFilteredRows,
    rows,
    setFilter,
    setHiddenColumns,
    state: { filters },
  } = tableInstance;
  const { categoryViews, onFilter, multiPanelSelectedUIState } =
    useCategoryFilter(
      preFilteredRows,
      categories,
      filters,
      setFilter,
      initialMultiPanelSelectedUIState
    );

  // Updates table hidden columns state.
  useEffect(() => {
    setHiddenColumns(hiddenColumns);
  }, [hiddenColumns, setHiddenColumns]);

  // Store latest filter state.
  useEffect(() => {
    storeFilters(filters);
    storeMultiPanelSelectedUIState(multiPanelSelectedUIState);
  }, [
    filters,
    multiPanelSelectedUIState,
    storeFilters,
    storeMultiPanelSelectedUIState,
  ]);

  // Handle side bar open/closed state beyond scope of component.
  const [isSideBarOpen, storeIsSideBarOpen] = useSessionStorage<boolean>(
    KEYS.SIDE_BAR_COLLECTIONS,
    true
  );

  return (
    <>
      <Head>
        <title>Collections - CZ CELLxGENE Discover</title>
      </Head>
      {isError ? null : shouldShowLoader ? (
        <Loader />
      ) : (
        shouldShowCollections && (
          <>
            <SideBar
              label="Filters"
              isOpen={isSideBarOpen}
              onToggle={storeIsSideBarOpen}
            >
              <CategoryFilters
                filters={partitionCategoryViews(categoryViews, mode)}
                onFilter={onFilter}
              />
            </SideBar>
            <View>
              {mode === VIEW_MODE.CURATOR && <CreateCollection />}
              {!rows || rows.length === 0 ? (
                <GridHero>
                  <h3>No Results</h3>
                  <p>There are no collections matching those filters.</p>
                </GridHero>
              ) : (
                <CollectionsGrid
                  mode={mode}
                  tableCountSummary={buildTableCountSummary(
                    rows,
                    preFilteredRows
                  )}
                  // @ts-expect-error -- revisit tableInstance typing
                  tableInstance={tableInstance}
                />
              )}
            </View>
            {/* May be added in the future after sign off */}
            {/* <BottomBanner /> */}
          </>
        )
      )}
    </>
  );
}

/**
 * Partitions the category views for collections and curator mode filtering.
 * @param categoryViews - View models of categories to display.
 * @param mode - View mode.
 * @returns partitioned category views.
 */
function partitionCategoryViews(
  categoryViews: CategoryView[],
  mode: VIEW_MODE
): CategoryView[][] {
  const partitionedValues: CategoryView[][] = [[]];
  for (const categoryView of categoryViews) {
    if (
      CATEGORY_FILTER_PARTITION_LIST[mode].includes(
        categoryView.categoryFilterId
      )
    ) {
      if (!partitionedValues[1]) {
        // Create the second partition if not already defined.
        partitionedValues.push([]);
      }
      // Assign category to the second partition.
      partitionedValues[1].push(categoryView);
    } else {
      // Assign category to the first partition.
      partitionedValues[0].push(categoryView);
    }
  }
  return partitionedValues;
}
