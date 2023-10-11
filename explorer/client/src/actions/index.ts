import { AnyAction } from "redux";
import type { Config } from "../globals";
import * as globals from "../globals";
import { AnnoMatrixLoader, AnnoMatrixObsCrossfilter } from "../annoMatrix";
import {
  catchErrorsWrap,
  doJsonRequest,
  dispatchNetworkErrorMessageToUser,
} from "../util/actionHelpers";
import { loadUserColorConfig } from "../util/stateManager/colorHelpers";
import { removeLargeDatasets } from "../util/stateManager/datasetMetadataHelpers";
import * as selnActions from "./selection";
import * as viewActions from "./viewStack";
import * as embActions from "./embedding";
import * as genesetActions from "./geneset";
import { AppDispatch, GetState } from "../reducers";
import { EmbeddingSchema, Field, Schema } from "../common/types/schema";
import { ConvertedUserColors } from "../reducers/colors";
import type { DatasetMetadata, Dataset, S3URI } from "../common/types/entities";
import { postExplainNewTab } from "../components/framework/toasters";
import { KEYS } from "../components/util/localStorage";
import {
  storageGetTransient,
  storageSetTransient,
} from "../components/util/transientLocalStorage";
import { selectIsUserStateDirty } from "../selectors/global";
import { DataframeValue, LabelArray, LabelIndex } from "../util/dataframe";
import { packDiffExPdu, DiffExMode, DiffExArguments } from "../util/diffexpdu";
import { track } from "../analytics";
import { EVENTS } from "../analytics/events";
import AnnoMatrix from "../annoMatrix/annoMatrix";

function setGlobalConfig(config: Config) {
  /**
   * Set any global run-time config not _exclusively_ managed by the config reducer.
   * This should only set fields defined in globals.globalConfig.
   */
  globals.globalConfig.maxCategoricalOptionsToDisplay =
    config?.parameters?.["max-category-items"] ??
    globals.globalConfig.maxCategoricalOptionsToDisplay;
}

/*
return promise fetching user-configured colors
*/
async function userColorsFetchAndLoad(
  dispatch: AppDispatch,
): Promise<{ type: string; userColors: ConvertedUserColors }> {
  return fetchJson<{ [category: string]: { [label: string]: string } }>(
    "colors",
  ).then((response) =>
    dispatch({
      type: "universe: user color load success",
      userColors: loadUserColorConfig(response),
    }),
  );
}

async function s3URIFetch(): Promise<S3URI> {
  return fetchJson<S3URI>("s3_uri");
}

async function schemaFetch(): Promise<{ schema: Schema }> {
  return fetchJson<{ schema: Schema }>("schema");
}

async function configFetchAndLoad(dispatch: AppDispatch): Promise<Config> {
  const response = await fetchJson<{ config: globals.Config }>("config");
  const config = { ...globals.configDefaults, ...response.config };

  setGlobalConfig(config);

  dispatch({
    type: "configuration load complete",
    config,
  });
  return config;
}

/**
 * Fetch dataset metadata and dispatch save to store, including portal URL returned in /config.
 * @param dispatch Function facilitating update of store.
 * @param oldPrefix API prefix with dataset path that dataset metadata lives on. (Not S3 URI)
 * @param config Response from config endpoint containing collection ID for the current dataset.
 */
async function datasetMetadataFetchAndLoad(
  dispatch: AppDispatch,
  oldPrefix: string,
  config: Config,
): Promise<void> {
  try {
    const datasetMetadataResponse = await fetchJson<{
      metadata: DatasetMetadata;
    }>("dataset-metadata", oldPrefix);

    // Create new dataset array with large datasets removed
    const { metadata: datasetMetadata } = datasetMetadataResponse;
    const datasets = removeLargeDatasets(
      datasetMetadata.collection_datasets,
      globals.DATASET_MAX_CELL_COUNT,
    );

    const { links } = config;
    dispatch({
      type: "dataset metadata load complete",
      datasetMetadata: {
        ...datasetMetadata,
        collection_datasets: datasets,
      },
      portalUrl: links["collections-home-page"],
    });
  } catch (error) {
    dispatch({
      type: "dataset metadata load error",
      error,
    });
  }
}

interface GeneInfoAPI {
  ncbi_url: string;
  name: string;
  synonyms: string[];
  summary: string;
  show_warning_banner: boolean;
}
/**
 * Fetch gene summary information
 * @param geneID ensembl ID corresponding to gene to search
 * @param gene human-readable name of gene
 */
async function fetchGeneInfo(
  geneID: DataframeValue,
  gene: string,
): Promise<GeneInfoAPI | undefined> {
  const response = await fetchJson<GeneInfoAPI>(
    `geneinfo?geneID=${geneID}&gene=${gene}`,
  );
  return response;
}

function prefetchEmbeddings(annoMatrix: AnnoMatrix) {
  /*
  prefetch requests for all embeddings
  */
  const { schema } = annoMatrix;
  const available = schema.layout.obs.map((v: EmbeddingSchema) => v.name);
  available.forEach((embName: EmbeddingSchema["name"]) =>
    annoMatrix.prefetch(Field.emb, embName, globals.numBinsEmb),
  );
}

/*
Application bootstrap
*/
const doInitialDataLoad = (): ((
  dispatch: AppDispatch,
  getState: GetState,
) => void) =>
  catchErrorsWrap(async (dispatch: AppDispatch) => {
    dispatch({ type: "initial data load start" });
    if (!globals.API) throw new Error("API not set");

    try {
      const s3URI = await s3URIFetch();
      const oldPrefix = globals.updateAPIWithS3(s3URI);
      const [config, schema] = await Promise.all([
        configFetchAndLoad(dispatch),
        schemaFetch(),
        userColorsFetchAndLoad(dispatch),
      ]);

      datasetMetadataFetchAndLoad(dispatch, oldPrefix, config);

      const baseDataUrl = `${globals.API.prefix}${globals.API.version}`;
      const annoMatrix = new AnnoMatrixLoader(baseDataUrl, schema.schema);
      const obsCrossfilter = new AnnoMatrixObsCrossfilter(annoMatrix);
      prefetchEmbeddings(annoMatrix);

      dispatch({
        type: "annoMatrix: init complete",
        annoMatrix,
        obsCrossfilter,
      });
      dispatch({ type: "initial data load complete" });

      const defaultEmbedding = config?.parameters?.default_embedding;
      const layoutSchema = schema?.schema?.layout?.obs ?? [];
      if (
        defaultEmbedding &&
        layoutSchema.some((s: EmbeddingSchema) => s.name === defaultEmbedding)
      ) {
        dispatch(embActions.layoutChoiceAction(defaultEmbedding));
      }
    } catch (error) {
      dispatch({ type: "initial data load error", error });
    }
  }, true);

function requestSingleGeneExpressionCountsForColoringPOST(
  gene: string,
): AnyAction {
  return {
    type: "color by expression",
    gene,
  };
}

const requestUserDefinedGene = (gene: string): AnyAction => ({
  type: "request user defined gene success",
  gene,
});

const dispatchDiffExpErrors = (
  dispatch: AppDispatch,
  response: Response,
): void => {
  switch (response.status) {
    case 403:
      dispatchNetworkErrorMessageToUser(
        "Too many cells selected for differential expression calculation - please make a smaller selection.",
      );
      break;
    case 501:
      dispatchNetworkErrorMessageToUser(
        "Differential expression is not implemented.",
      );
      break;
    default: {
      const msg = `Unexpected differential expression HTTP response ${response.status}, ${response.statusText}`;
      dispatchNetworkErrorMessageToUser(msg);
      dispatch({
        type: "request differential expression error",
        error: new Error(msg),
      });
    }
  }
};

const requestDifferentialExpression =
  (set1: LabelArray, set2: LabelArray, num_genes = 15) =>
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types -- TODO: type diff exp data in genesets reducer (client/src/reducers/genesets.ts)
  async (dispatch: AppDispatch, getState: GetState) => {
    // (thuang): For measuring performance
    const startMs = performance.now();
    const cellCount1 = set1.length;
    const cellCount2 = set2.length;

    dispatch({ type: "request differential expression started" });
    try {
      if (!globals.API) throw new Error("API not set");
      /*
      Steps:
      1. get the most differentially expressed genes
      2. get expression data for each
      */
      const { annoMatrix } = getState();
      const varIndexName = annoMatrix.schema.annotations.var.index;

      // // Legal values are null, Array or TypedArray.  Null is initial state.
      if (!set1) set1 = new Int32Array();
      if (!set2) set2 = new Int32Array();

      // set1/set2 are LabelArray, which may be Array<> or Int32Array.
      // The API accepts Uint32Array.
      const set1Uint32 = new Uint32Array(set1 as Int32Array);
      const set2Uint32 = new Uint32Array(set2 as Int32Array);
      const deArgs: DiffExArguments = {
        mode: DiffExMode.TopN,
        params: { N: num_genes },
        set1: set1Uint32,
        set2: set2Uint32,
      };

      const res = await fetch(
        `${globals.API.prefix}${globals.API.version}diffexp/obs2`,
        {
          method: "POST",
          headers: new Headers({
            Accept: "application/json",
            "Content-Type": "application/octet-stream",
          }),
          body: packDiffExPdu(deArgs),
          credentials: "include",
        },
      );

      if (!res.ok || res.headers.get("Content-Type") !== "application/json") {
        track(EVENTS.EXPLORER_DIFF_EXP_BUTTON_CLICKED, {
          cellCount1,
          cellCount2,
          timeToComplete: performance.now() - startMs,
          status: "error",
        });

        return dispatchDiffExpErrors(dispatch, res);
      }

      const response = await res.json();
      const varIndex = await annoMatrix.fetch("var", varIndexName);
      const diffexpLists = { negative: [], positive: [] };
      for (const polarity of Object.keys(
        diffexpLists,
      ) as (keyof typeof diffexpLists)[]) {
        diffexpLists[polarity] = response[polarity].map(
          // TODO: swap out with type defined at genesets reducer when made
          (v: [LabelIndex, number, number, number]) => [
            varIndex.at(v[0], varIndexName),
            ...v.slice(1),
          ],
        );
      }

      track(EVENTS.EXPLORER_DIFF_EXP_BUTTON_CLICKED, {
        cellCount1,
        cellCount2,
        timeToComplete: performance.now() - startMs,
        status: "success",
      });

      /* then send the success case action through */
      return dispatch({
        type: "request differential expression success",
        data: diffexpLists,
      });
    } catch (error) {
      track(EVENTS.EXPLORER_DIFF_EXP_BUTTON_CLICKED, {
        cellCount1,
        cellCount2,
        timeToComplete: performance.now() - startMs,
        status: "error",
      });

      return dispatch({
        type: "request differential expression error",
        error,
      });
    }
  };

/**
 * Check local storage for flag indicating that the work in progress toast should be displayed.
 */
const checkExplainNewTab =
  () =>
  (dispatch: AppDispatch): void => {
    const workInProgressWarn = storageGetTransient(KEYS.WORK_IN_PROGRESS_WARN);
    if (workInProgressWarn) {
      dispatch({ type: "work in progress warning displayed" });
      postExplainNewTab(
        "To maintain your in-progress work on the previous dataset, we opened this dataset in a new tab.",
      );
    }
  };

/**
 * Navigate to URL in the same browser tab if there is no work in progress, otherwise open URL in new tab.
 * @param url - URL to navigate to.
 */
const navigateCheckUserState =
  (url: string): ((dispatch: AppDispatch, getState: GetState) => void) =>
  (_dispatch: AppDispatch, getState: GetState) => {
    const workInProgress = selectIsUserStateDirty(getState());
    if (workInProgress) {
      openTab(`${url}?${globals.QUERY_PARAM_EXPLAIN_NEW_TAB}`);
    } else {
      window.location.href = url;
    }
  };

/**
 * Handle select of dataset from dataset selector: determine whether to display dataset in current browser tab or open
 * dataset in new tab if user currently has work in progress.
 * @param dataset Dataset to switch to and load in the current tab.
 */
const selectDataset =
  (dataset: Dataset): ((dispatch: AppDispatch, getState: GetState) => void) =>
  (dispatch: AppDispatch, getState: GetState) => {
    const workInProgress = selectIsUserStateDirty(getState());
    if (workInProgress) {
      dispatch(openDataset(dataset));
    } else {
      dispatch(switchDataset(dataset));
    }
  };

/**
 * Open selected dataset in a new tab. Create local storage with expiry to pop toast once dataset is opened.
 * @param dataset Dataset to open in new tab.
 */
const openDataset =
  (dataset: Dataset): ((dispatch: AppDispatch) => void) =>
  (dispatch: AppDispatch) => {
    const deploymentUrl = dataset.dataset_deployments?.[0].url;
    if (!deploymentUrl) {
      dispatchNetworkErrorMessageToUser("Unable to open dataset.");
      return;
    }

    dispatch({ type: "dataset opened" });
    storageSetTransient(KEYS.WORK_IN_PROGRESS_WARN, 10000);
    openTab(deploymentUrl);
  };

/**
 * Open new tab and navigate to the given URL.
 * @param url - URL to navigate to.
 */
const openTab = (url: string) => {
  window.open(url, "_blank", "noopener");
};

/**
 * Open selected dataset in a new tab.
 * @param dataset Dataset to open in new tab.
 */
const switchDataset =
  (dataset: Dataset): ((dispatch: AppDispatch) => void) =>
  (dispatch: AppDispatch) => {
    dispatch({ type: "reset" });
    dispatch({ type: "dataset switch" });

    const deploymentUrl = dataset.dataset_deployments?.[0].url;
    if (!deploymentUrl) {
      dispatchNetworkErrorMessageToUser("Unable to switch datasets.");
      return;
    }
    dispatch(updateLocation(deploymentUrl));
    globals.updateApiPrefix();
    dispatch(doInitialDataLoad());
  };

/**
 * Update browser location by adding corresponding entry to the session's history stack.
 * @param url - URL to update browser location to.
 */
const updateLocation = (url: string) => (dispatch: AppDispatch) => {
  dispatch({ type: "location update" });
  window.history.pushState(null, "", url);
};

function fetchJson<T>(pathAndQuery: string, apiPrefix?: string): Promise<T> {
  if (!globals.API) throw new Error("API not initialized");
  if (!apiPrefix) apiPrefix = globals.API.prefix;
  return doJsonRequest<T>(
    `${apiPrefix}${globals.API.version}${pathAndQuery}`,
  ) as Promise<T>;
}

export default {
  doInitialDataLoad,
  requestDifferentialExpression,
  requestSingleGeneExpressionCountsForColoringPOST,
  requestUserDefinedGene,
  checkExplainNewTab,
  navigateCheckUserState,
  selectDataset,
  fetchGeneInfo,
  selectContinuousMetadataAction: selnActions.selectContinuousMetadataAction,
  selectCategoricalMetadataAction: selnActions.selectCategoricalMetadataAction,
  selectCategoricalAllMetadataAction:
    selnActions.selectCategoricalAllMetadataAction,
  graphBrushStartAction: selnActions.graphBrushStartAction,
  graphBrushChangeAction: selnActions.graphBrushChangeAction,
  graphBrushDeselectAction: selnActions.graphBrushDeselectAction,
  graphBrushCancelAction: selnActions.graphBrushCancelAction,
  graphBrushEndAction: selnActions.graphBrushEndAction,
  graphLassoStartAction: selnActions.graphLassoStartAction,
  graphLassoEndAction: selnActions.graphLassoEndAction,
  graphLassoCancelAction: selnActions.graphLassoCancelAction,
  graphLassoDeselectAction: selnActions.graphLassoDeselectAction,
  clipAction: viewActions.clipAction,
  subsetAction: viewActions.subsetAction,
  resetSubsetAction: viewActions.resetSubsetAction,
  layoutChoiceAction: embActions.layoutChoiceAction,
  setCellSetFromSelection: selnActions.setCellSetFromSelection,
  genesetDelete: genesetActions.genesetDelete,
  genesetAddGenes: genesetActions.genesetAddGenes,
  genesetDeleteGenes: genesetActions.genesetDeleteGenes,
};
