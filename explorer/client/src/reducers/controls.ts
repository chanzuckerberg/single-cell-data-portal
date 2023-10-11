import { AnyAction } from "redux";

type Level = "top" | "bottom" | "";

interface StackLevels {
  geneLevel: Level;
  scatterplotLevel: Level;
}
/* logic for minimizing and maximizing pop-ups */
const minimizeMaximizePopUps = (
  geneLevel: Level,
  geneIsMinimized: boolean,
  geneIsOpen: boolean,
  scatterplotLevel: Level,
  scatterplotIsMinimized: boolean,
  scatterplotIsOpen: boolean,
): StackLevels => {
  if (
    geneIsMinimized &&
    geneIsOpen &&
    scatterplotIsMinimized &&
    scatterplotIsOpen
  ) {
    return {
      geneLevel,
      scatterplotLevel,
    };
  }
  if (!geneIsMinimized) {
    return {
      geneLevel: scatterplotIsMinimized && scatterplotIsOpen ? "top" : "bottom",
      scatterplotLevel:
        scatterplotIsMinimized && scatterplotIsOpen ? "bottom" : "",
    };
  }
  if (!scatterplotIsMinimized) {
    return {
      scatterplotLevel: geneIsMinimized && geneIsOpen ? "top" : "bottom",
      geneLevel: geneIsMinimized && geneIsOpen ? "bottom" : "",
    };
  }
  return {
    geneLevel: geneIsMinimized ? "bottom" : "top",
    scatterplotLevel: scatterplotIsMinimized ? "bottom" : "top",
  };
};

interface ControlsState {
  loading: boolean;
  error: Error | string | null;
  resettingInterface: boolean;
  graphInteractionMode: "zoom" | "select";
  opacityForDeselectedCells: number;
  scatterplotXXaccessor: string | false;
  scatterplotYYaccessor: string | false;
  geneIsOpen: boolean;
  scatterplotIsMinimized: boolean;
  geneIsMinimized: boolean;
  scatterplotLevel: Level;
  scatterplotIsOpen: boolean;
  geneLevel: Level;
  gene: string | null;
  infoError: string | null;
  graphRenderCounter: number;
  colorLoading: boolean;
  datasetDrawer: boolean;
  geneUrl: string;
  geneSummary: string;
  geneName: string;
  geneSynonyms: string[];
}
const Controls = (
  state: ControlsState = {
    // data loading flag
    loading: true,
    error: null,
    // all of the data + selection state
    resettingInterface: false,
    graphInteractionMode: "select",
    opacityForDeselectedCells: 0.2,
    scatterplotXXaccessor: false, // just easier to read
    scatterplotYYaccessor: false,
    geneIsOpen: false,
    scatterplotIsMinimized: false,
    geneIsMinimized: false,
    geneUrl: "",
    geneSummary: "",
    geneSynonyms: [""],
    geneName: "",
    scatterplotLevel: "",
    scatterplotIsOpen: false,
    geneLevel: "",
    gene: null,
    infoError: null,
    graphRenderCounter: 0 /* integer as <Component key={graphRenderCounter} - a change in key forces a remount */,
    colorLoading: false,
    datasetDrawer: false,
  },
  action: AnyAction,
) => {
  /*
  For now, log anything looking like an error to the console.
  */

  if (action.error || /error/i.test(action.type)) {
    console.error(action.error);
  }

  switch (action.type) {
    case "initial data load start": {
      return { ...state, loading: true };
    }
    case "initial data load complete": {
      /* now fully loaded */
      return {
        ...state,
        loading: false,
        error: null,
        resettingInterface: false,
      };
    }
    case "reset subset": {
      return {
        ...state,
        resettingInterface: false,
      };
    }
    case "subset to selection": {
      return {
        ...state,
        loading: false,
        error: null,
      };
    }
    case "initial data load error": {
      return {
        ...state,
        loading: false,
        error: action.error,
      };
    }
    case "color by geneset mean expression": {
      return {
        ...state,
        colorLoading: true,
      };
    }
    case "color by geneset mean expression success": {
      return {
        ...state,
        colorLoading: false,
      };
    }
    /*******************************
             User Events
     *******************************/
    case "change graph interaction mode":
      return {
        ...state,
        graphInteractionMode: action.data,
      };
    case "change opacity deselected cells in 2d graph background":
      return {
        ...state,
        opacityForDeselectedCells: action.data,
      };
    case "increment graph render counter": {
      const c = state.graphRenderCounter + 1;
      return {
        ...state,
        graphRenderCounter: c,
      };
    }

    /*******************************
              Gene Info
    *******************************/
    case "open gene info": {
      // if scatterplot is already open
      if (
        !state.scatterplotIsMinimized &&
        state.scatterplotXXaccessor &&
        state.scatterplotYYaccessor
      ) {
        state.scatterplotIsMinimized = true;
        state.geneIsMinimized = false;
      }

      const stackLevels = minimizeMaximizePopUps(
        state.geneLevel,
        state.geneIsMinimized,
        true,
        state.scatterplotLevel,
        state.scatterplotIsMinimized,
        state.scatterplotIsOpen,
      );

      return {
        ...state,
        geneIsOpen: true,
        gene: action.gene,
        geneUrl: action.url,
        geneSummary: action.summary,
        geneSynonyms: action.synonyms,
        geneName: action.name,
        geneIsMinimized: false,
        geneLevel: stackLevels.geneLevel,
        infoError: action.infoError,
        showWarningBanner: action.showWarningBanner,
        scatterplotLevel: stackLevels.scatterplotLevel,
      };
    }

    case "load gene info": {
      // if scatterplot is already open
      if (
        !state.scatterplotIsMinimized &&
        state.scatterplotXXaccessor &&
        state.scatterplotYYaccessor
      ) {
        state.scatterplotIsMinimized = true;
        state.geneIsMinimized = false;
      }

      const stackLevels = minimizeMaximizePopUps(
        state.geneLevel,
        state.geneIsMinimized,
        true,
        state.scatterplotLevel,
        state.scatterplotIsMinimized,
        state.scatterplotIsOpen,
      );

      return {
        ...state,
        geneIsOpen: true,
        gene: action.gene,
        geneUrl: "",
        geneSummary: "",
        geneSynonyms: [""],
        geneName: "",
        geneIsMinimized: false,
        geneLevel: stackLevels.geneLevel,
        infoError: null,
        scatterplotLevel: stackLevels.scatterplotLevel,
      };
    }

    case "minimize/maximize gene info": {
      state.geneIsMinimized = !state.geneIsMinimized;
      if (!state.geneIsMinimized) {
        state.scatterplotIsMinimized = true;
      }
      const stackLevels = minimizeMaximizePopUps(
        state.geneLevel,
        state.geneIsMinimized,
        state.geneIsOpen,
        state.scatterplotLevel,
        state.scatterplotIsMinimized,
        state.scatterplotIsOpen,
      );

      return {
        ...state,
        geneIsMinimized: state.geneIsMinimized,
        geneLevel: stackLevels.geneLevel,
        scatterplotIsMinimized: state.scatterplotIsMinimized,
        scatterplotLevel: stackLevels.scatterplotLevel,
        infoError: state.infoError,
      };
    }

    case "clear gene info": {
      const stackLevels = minimizeMaximizePopUps(
        state.geneLevel,
        state.geneIsMinimized,
        false,
        state.scatterplotLevel,
        state.scatterplotIsMinimized,
        state.scatterplotIsOpen,
      );

      return {
        ...state,
        geneIsOpen: false,
        geneIsMinimized: null,
        geneLevel: stackLevels.geneLevel,
        scatterplotLevel: stackLevels.scatterplotLevel,
        infoError: action.infoError,
      };
    }

    /*******************************
              Scatterplot
    *******************************/
    case "set scatterplot x": {
      // gene info is already open
      if (
        !state.geneIsMinimized &&
        state.geneIsOpen &&
        state.scatterplotYYaccessor
      ) {
        state.geneIsMinimized = true;
        state.scatterplotIsMinimized = false;
        state.scatterplotIsOpen = true;
      }

      const stackLevels = minimizeMaximizePopUps(
        state.geneLevel,
        state.geneIsMinimized,
        state.geneIsOpen,
        state.scatterplotLevel,
        state.scatterplotIsMinimized,
        state.scatterplotIsOpen,
      );
      state.geneLevel = stackLevels.geneLevel;
      state.scatterplotLevel = stackLevels.scatterplotLevel;

      return {
        ...state,
        scatterplotXXaccessor: action.data,
        scatterplotIsMinimized: false,
        scatterplotLevel: state.scatterplotLevel,
      };
    }

    case "set scatterplot y": {
      // gene info is already open
      if (
        !state.geneIsMinimized &&
        state.geneIsOpen &&
        state.scatterplotXXaccessor
      ) {
        state.geneIsMinimized = true;
        state.scatterplotIsMinimized = false;
        state.scatterplotIsOpen = true;
      }

      const stackLevels = minimizeMaximizePopUps(
        state.geneLevel,
        state.geneIsMinimized,
        state.geneIsOpen,
        state.scatterplotLevel,
        state.scatterplotIsMinimized,
        state.scatterplotIsOpen,
      );
      state.geneLevel = stackLevels.geneLevel;
      state.scatterplotLevel = stackLevels.scatterplotLevel;

      return {
        ...state,
        scatterplotYYaccessor: action.data,
        scatterplotIsMinimized: false,
        scatterplotLevel: state.scatterplotLevel,
      };
    }

    case "minimize/maximize scatterplot": {
      state.scatterplotIsMinimized = !state.scatterplotIsMinimized;
      if (!state.scatterplotIsMinimized) {
        state.geneIsMinimized = true;
      }
      const stackLevels = minimizeMaximizePopUps(
        state.geneLevel,
        state.geneIsMinimized,
        state.geneIsOpen,
        state.scatterplotLevel,
        state.scatterplotIsMinimized,
        state.scatterplotIsOpen,
      );
      state.geneLevel = stackLevels.geneLevel;
      state.scatterplotLevel = stackLevels.scatterplotLevel;

      return {
        ...state,
        scatterplotIsMinimized: state.scatterplotIsMinimized,
        scatterplotLevel: state.scatterplotLevel,
        geneIsMinimized: state.geneIsMinimized,
      };
    }

    case "clear scatterplot": {
      state.scatterplotIsOpen = false;
      const stackLevels = minimizeMaximizePopUps(
        state.geneLevel,
        state.geneIsMinimized,
        state.geneIsOpen,
        state.scatterplotLevel,
        state.scatterplotIsMinimized,
        state.scatterplotIsOpen,
      );
      state.geneLevel = stackLevels.geneLevel;
      state.scatterplotLevel = stackLevels.scatterplotLevel;

      return {
        ...state,
        scatterplotXXaccessor: null,
        scatterplotYYaccessor: null,
        scatterplotIsMinimized: null,
        scatterplotLevel: state.scatterplotLevel,
      };
    }

    /**************************
          Dataset Drawer
     **************************/
    case "toggle dataset drawer":
      return { ...state, datasetDrawer: !state.datasetDrawer };

    default:
      return state;
  }
};

export default Controls;
