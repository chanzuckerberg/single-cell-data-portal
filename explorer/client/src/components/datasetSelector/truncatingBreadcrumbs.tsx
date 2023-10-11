/* Core dependencies */
import { BreadcrumbProps, Classes, ResizeSensor } from "@blueprintjs/core";
import { ResizeEntry } from "@blueprintjs/core/src/components/resize-sensor/resizeObserverTypes";
import React, { CSSProperties, useCallback, useEffect, useState } from "react";

/* App dependencies */
import Truncate from "../util/truncate";

interface Props {
  breadcrumbRenderer: (item: TruncatingBreadcrumbProps) => JSX.Element;
  currentBreadcrumbRenderer: (item: TruncatingBreadcrumbProps) => JSX.Element;
  items: TruncatingBreadcrumbProps[];
}

/**
 * Internal representation of breadcrumb props, with additional config to help calculate truncated state, if any.
 */
interface BreadcrumbConfig {
  displayAs: React.ReactNode;
  hidden: boolean;
  item: TruncatingBreadcrumbProps;
  fullWidth: number; // Width of breadcrumb when in F (full text) state
  shortTextWidth: number; // Width of breadcrumb when in S (short text) state
  truncatedWidth: number; // Width of breadcrumb when in T (truncated text) state
}

/**
 * Function used to determine actual rendered widths of text.
 */
type MeasureFn = (text: string) => number;

/**
 * Standard breadcrumb link.
 */
export interface TruncatingBreadcrumbProps extends BreadcrumbProps {
  shortText: string;
}

// Characters to be used to indicate display text has been truncated
const CHAR_ELLIPSIS = "...";

// Width of breadcrumb icon, if specified.
const ICON_WIDTH = 26;

// Approximate padding in pixels for each breadcrumb, including caret.
const ITEM_PADDING = 26;

// Minimum number of characters to be displayed before transitioning to a smaller state of the breadcrumbs
const MIN_VISIBLE_CHARS = 11;

// Inline styles for breadcrumbs list element
const STYLE_BREADCRUMBS: CSSProperties = {
  display: "flex",
  flexWrap: "nowrap",
  whiteSpace: "nowrap",
};

// Inline styles for hidden measure element
const STYLE_MEASURE_SPAN: CSSProperties = {
  left: "-999px",
  position: "absolute",
  top: "-999px",
  visibility: "hidden",
};

/**
 * Individual Breadcrumb States
 * ----------------------------
 * F - full text
 * T - truncated text
 * S - indicates use of short text (eg "Collection" for collection name or "Dataset" for dataset name)
 * H - hidden
 */
enum BreadcrumbState {
  "FULL" = "F", // eg "Tabula Muris Senis"
  "TRUNCATED" = "T", // eg "Tabula...Senis"
  "SHORT_TEXT" = "S", // eg "Collection"
  "HIDDEN" = "H", // --
}

/**
 * Breadcrumbs State
 * -----------------
 * FFF - full, full, full
 * FTF - full, truncated, full
 * HSF - hidden, short text, full
 * HST - hidden, short text, truncated
 * HHS - hidden, hidden, short text
 * HHH - hidden, hidden, hidden
 */
type BreadcrumbsState = [BreadcrumbState, BreadcrumbState, BreadcrumbState];
const STATES_FFF: BreadcrumbsState = [
  BreadcrumbState.FULL,
  BreadcrumbState.FULL,
  BreadcrumbState.FULL,
];
const STATES_FTF: BreadcrumbsState = [
  BreadcrumbState.FULL,
  BreadcrumbState.TRUNCATED,
  BreadcrumbState.FULL,
];
const STATES_HSF: BreadcrumbsState = [
  BreadcrumbState.HIDDEN,
  BreadcrumbState.SHORT_TEXT,
  BreadcrumbState.FULL,
];
const STATES_HST: BreadcrumbsState = [
  BreadcrumbState.HIDDEN,
  BreadcrumbState.SHORT_TEXT,
  BreadcrumbState.TRUNCATED,
];
const STATES_HHS: BreadcrumbsState = [
  BreadcrumbState.HIDDEN,
  BreadcrumbState.HIDDEN,
  BreadcrumbState.SHORT_TEXT,
];
const STATES_HHH: BreadcrumbsState = [
  BreadcrumbState.HIDDEN,
  BreadcrumbState.HIDDEN,
  BreadcrumbState.HIDDEN,
];

/**
 * Breadcrumb Transitions
 * ----------------------
 * Breadcrumbs can transition bidirectionally through states in the following order, and can also repeat individual states.
 */
const BreadcrumbsStateTransitions: BreadcrumbsState[] = [
  STATES_FFF,
  STATES_FTF,
  STATES_HSF,
  STATES_HST,
  STATES_HHS,
  STATES_HHH,
];

/**
 * Tuple of links rendered with cellxgene-specific truncating functionality.
 */
const TruncatingBreadcrumbs = React.memo<Props>(
  ({ breadcrumbRenderer, currentBreadcrumbRenderer, items }) => {
    // Total width available to truncating breadcrumbs, set from ResizeSensor's resize event.
    const [availableWidth, setAvailableWidth] = useState<number>(0);

    // View models backing calculations and rendering of truncating breadcrumbs.
    const [breadcrumbConfigs, setBreadcrumbConfigs] = useState<
      BreadcrumbConfig[]
    >([]);

    // Flag indicating fonts have been loaded and calculations of widths for each breadcrumb state for each breadcrumb
    // can begin.
    const [fontLoaded, setFontLoaded] = useState<boolean>(false);

    // "Measuring tape" element; used to calculate width of breadcrumb text in each breadcrumb state. Total number of
    // calculations is n * 3 where n is the number of breadcrumbs.
    const [measureEl, setMeasureEl] = useState<HTMLSpanElement>();

    // Set measure element once rendered.
    const measuredRef = useCallback((el: HTMLSpanElement) => {
      if (el !== null) {
        setMeasureEl(el);
      }
    }, []);

    // Load of breadcrumbs can begin once fonts are loaded.
    useEffect(() => {
      document.fonts.ready.then(() => {
        setFontLoaded(true);
      });
    }, []);

    // Once we have the available width, measure element and font, build up view model.
    useEffect(() => {
      if (availableWidth > 0 && fontLoaded && items.length > 0 && measureEl) {
        setBreadcrumbConfigs(
          buildBreadcrumbConfigs(
            items,
            availableWidth,
            measureTextWidth(measureEl),
          ),
        );
      }
    }, [availableWidth, fontLoaded, items, measureEl]);

    /**
     * On resize callback from ResizeSensor, save the current width of the breadcrumbs.
     * @param entries - Array of elements being observed for resize events.
     */
    const onResize = (entries: ResizeEntry[]) => {
      setAvailableWidth(Math.floor(entries[0].contentRect.width));
    };

    /**
     * Build list elements for each breadcrumb config.
     * @param configs - The set of config objects backing each breadcrumb.
     * @returns List elements for each breadcrumb config.
     */
    const renderBreadcrumbs = (
      configs: BreadcrumbConfig[],
    ): (JSX.Element | null)[] =>
      configs.map((config: BreadcrumbConfig, i: number) => {
        if (config.hidden) {
          return null;
        }
        // Create new version of breadcrumb props from original item with updated text for display.
        const props = {
          ...config.item,
          text: config.displayAs,
        };
        return (
          <li
            data-testid={`bc-${config.item.shortText}`}
            key={`bc-${config.item.shortText}`}
          >
            {isCurrentBreadcrumb(configs, i)
              ? currentBreadcrumbRenderer(props)
              : breadcrumbRenderer(props)}
          </li>
        );
      });

    return (
      <div
        style={{
          overflow: "visible",
          position: "absolute",
          bottom: 8,
          left: 8,
          width: "100%",
        }}
      >
        <ResizeSensor onResize={onResize}>
          <ul className={Classes.BREADCRUMBS} style={STYLE_BREADCRUMBS}>
            {renderBreadcrumbs(breadcrumbConfigs)}
          </ul>
        </ResizeSensor>
        <span ref={measuredRef} style={STYLE_MEASURE_SPAN} />
      </div>
    );
  },
);

/**
 * Build breadcrumb view models to fit the given available width.
 * @param items - Set of breadcrumb props passed in via props.
 * @param availableWidth - Full width available to breadcrumbs.
 * @param measureFn - Function called to measure actual rendered pixel width of text.
 * @returns The state for each breadcrumb given the width available. For example, [F, F, F].
 */
function buildBreadcrumbConfigs(
  items: TruncatingBreadcrumbProps[],
  availableWidth: number,
  measureFn: MeasureFn,
): BreadcrumbConfig[] {
  // Build default config object for each breadcrumb including measured widths.
  const configs = buildDefaultBreadcrumbConfigs(items, measureFn);

  // Calculate the state of each breadcrumb for the available width. For example, [F, F, F].
  const breadcrumbsState = calculateBreadcrumbsStateForAvailableWidth(
    configs,
    availableWidth,
  );

  // Update breadcrumb configs to match the calculated breadcrumbs state for the available width.
  return configs.map((configDefaults: BreadcrumbConfig, i: number) => {
    const breadcrumbState = breadcrumbsState[i];
    const displayAs = getBreadcrumbDisplayAs(
      configs,
      i,
      breadcrumbsState,
      availableWidth,
    );
    return {
      ...configDefaults,
      displayAs,
      hidden: isStateHidden(breadcrumbState),
    };
  });
}

/**
 * Build initial view model for each breadcrumb.
 * @param items - Set of breadcrumb props passed in via props.
 * @param measureFn - Function called to measure actual rendered pixel width of text.
 * @returns Breadcrumb configs set with actual DOM size calculations and defaults for all configuration values.
 */
function buildDefaultBreadcrumbConfigs(
  items: TruncatingBreadcrumbProps[],
  measureFn: MeasureFn,
): BreadcrumbConfig[] {
  return items.map((item: TruncatingBreadcrumbProps, i: number) => {
    // Determine any additional width required by the breadcrumb, either width for an icon or caret, or both.
    const iconWidth = item.icon ? ICON_WIDTH : 0;
    const padding = isCurrentBreadcrumb(items, i) ? 0 : ITEM_PADDING;
    const baseWidth = iconWidth + padding;

    // Determine minimum visible text in truncated state.
    const text = item.text as string;
    const truncatedText = truncateToLength(text, MIN_VISIBLE_CHARS);

    return {
      displayAs: text,
      fullWidth: measureFn(text) + baseWidth,
      hidden: false,
      item,
      shortTextWidth: measureFn(item.shortText) + baseWidth,
      truncatedWidth: measureFn(truncatedText) + baseWidth,
    };
  });
}

/**
 * Return the width that the truncated item has available for display. That is, the available width minus the widths
 * required by the other, non-truncated, items.
 * @param configs - Set of config objects backing each breadcrumb.
 * @param truncatedIndex - Index of breadcrumb currently being truncated.
 * @param breadcrumbsState - State of each breadcrumb. For example, [F, F, F].
 * @param availableWidth - Full width available to breadcrumbs.
 * @returns Width that truncated breadcrumb has available for display.
 */
function calculateAvailableTruncatedWidth(
  configs: BreadcrumbConfig[],
  truncatedIndex: number,
  breadcrumbsState: BreadcrumbsState,
  availableWidth: number,
): number {
  // Grab the configs other than the truncated config.
  const otherItems = [...configs];
  otherItems.splice(truncatedIndex, 1);

  // Grab the states of the configs other than the truncated config.
  const otherItemsState = [...breadcrumbsState];
  otherItemsState.splice(truncatedIndex, 1);

  // Calculate the width of the other items in their corresponding states.
  const otherItemsRequiredWidth = calculateRequiredWidth(
    otherItemsState,
    otherItems,
  );

  // Subtract icon width of breadcrumb, if icon specified.
  const { item } = configs[truncatedIndex];
  const iconWidth = item.icon ? ICON_WIDTH : 0;

  // Subtract padding if icon is specified
  const paddingWidth = item.icon ? 0 : ITEM_PADDING;

  return availableWidth - otherItemsRequiredWidth - iconWidth - paddingWidth;
}

/**
 * Determine the current items state (eg FFF, FTF etc) for the given available width and set of items.
 * @param configs - Set of config objects backing each breadcrumb.
 * @param availableWidth - Full width available to breadcrumbs.
 * @returns The ideal state for each breadcrumb given the width available. For example, [F, F, F].
 */
function calculateBreadcrumbsStateForAvailableWidth(
  configs: BreadcrumbConfig[],
  availableWidth: number,
): BreadcrumbsState {
  for (let i = 0; i < BreadcrumbsStateTransitions.length; i += 1) {
    const itemsState = BreadcrumbsStateTransitions[i];
    const requiredWidth = calculateRequiredWidth(itemsState, configs);
    if (availableWidth >= requiredWidth) {
      return itemsState;
    }
  }
  return BreadcrumbsStateTransitions[BreadcrumbsStateTransitions.length - 1]; // There's a problem, default to smallest state.
}

/**
 * Calculate the total width required to display the given breadcrumbs in the given breadcrumb states.
 * @param breadcrumbsState - State of each breadcrumb that width required to display is being calculated for. For
 * example, [F, F, F].
 * @param configs - Set of config objects backing each breadcrumb.
 * @returns Number presenting the total number of pixels required to display the items in their current states.
 */
function calculateRequiredWidth(
  breadcrumbsState: BreadcrumbState[],
  configs: BreadcrumbConfig[],
): number {
  return configs.reduce(
    (accum: number, config: BreadcrumbConfig, i: number) => {
      // Grab the state for this breadcrumb. For example, given the state HTF, the state of the first item is H, the state
      // of the second item is T and the state of the third item is F.
      const itemState = breadcrumbsState[i];
      if (!itemState) {
        return accum;
      }
      // Add the width of text corresponding to the item's state.
      if (isStateShortText(itemState)) {
        accum += config.shortTextWidth;
      } else if (isStateTruncated(itemState)) {
        accum += config.truncatedWidth;
      } else if (isStateFull(itemState)) {
        accum += config.fullWidth;
      }
      return accum;
    },
    0,
  );
}

/**
 * Return the display for the given breadcrumb in the given state.
 * @returns Short text or full text if breadcrumb is in corresponding short text/full text mode, otherwise truncate component
 */
function getBreadcrumbDisplayAs(
  configs: BreadcrumbConfig[],
  breadcrumbIndex: number,
  breadcrumbsState: BreadcrumbsState,
  availableWidth: number,
): React.ReactNode {
  const breadcrumbState = breadcrumbsState[breadcrumbIndex];
  const config = configs[breadcrumbIndex];
  if (isStateShortText(breadcrumbState)) {
    return config.item.shortText;
  }
  if (isStateTruncated(breadcrumbState)) {
    const truncatedAvailableWidth = calculateAvailableTruncatedWidth(
      configs,
      breadcrumbIndex,
      breadcrumbsState,
      availableWidth,
    );
    return (
      <Truncate>
        <span
          style={{
            display: "block",
            width: `${truncatedAvailableWidth}px`,
          }}
        >
          {config.item.text}
        </span>
      </Truncate>
    );
  }
  return config.item.text as string;
}

/**
 * Determine if breadcrumb is the current breadcrumb.
 * @param items - Set of objects backing each breadcrumb.
 * @param index - Position of breadcrumb in set.
 * @returns True if breadcrumb is the last breadcrumb in the set of breadcrumbs.
 */
function isCurrentBreadcrumb(
  items: BreadcrumbConfig[] | TruncatingBreadcrumbProps[],
  index: number,
): boolean {
  return index === items.length - 1;
}

/**
 * Returns true if state is full.
 * @param state - Name of state to check.
 * @returns True if given state is full.
 */
function isStateFull(state: BreadcrumbState): boolean {
  return state === BreadcrumbState.FULL;
}

/**
 * Returns true if state is hidden.
 * @param state - Name of state to check.
 * @returns True if given state is hidden.
 */
function isStateHidden(state: BreadcrumbState): boolean {
  return state === BreadcrumbState.HIDDEN;
}

/**
 * Returns true if state is short text.
 * @param state - Name of state to check.
 * @returns True if given state is short text.
 */
function isStateShortText(state: BreadcrumbState): boolean {
  return state === BreadcrumbState.SHORT_TEXT;
}

/**
 * Returns true if state is truncated.
 * @param state - Name of state to check.
 * @returns True if given state is truncated.
 */
function isStateTruncated(state: BreadcrumbState): boolean {
  return state === BreadcrumbState.TRUNCATED;
}

/**
 * Create function that returns actual rendered pixel width of given text.
 * @param measureEl - Rendered, hidden, element to update and measure.
 * @returns Function that returns the rendered pixel width of its given text.
 */
function measureTextWidth(measureEl: HTMLSpanElement): MeasureFn {
  return (text: string) => {
    measureEl.innerText = String(text);
    return measureEl.offsetWidth;
  };
}

/**
 * Return text truncated from center with ellipsis added.
 * @param text - Text to truncate
 * @param length - Length to truncate text to.
 * @returns Truncated text in format "start...end".
 */
function truncateToLength(text: string, length: number): string {
  if (text.length <= length) {
    return text;
  }
  // Determine the break indices for the "before" and "after" ellipsis text tokens
  const tokenBeforeEndIndex = Math.ceil(length / 2);
  const tokenAfterStartIndex = Math.floor(length / 2);
  // Split text at break indices and join with ellipsis
  const tokenBefore = text.substr(0, tokenBeforeEndIndex).trim();
  const tokenAfter = text.substr(text.length - tokenAfterStartIndex).trim();
  return `${tokenBefore}${CHAR_ELLIPSIS}${tokenAfter}`;
}

export default TruncatingBreadcrumbs;
