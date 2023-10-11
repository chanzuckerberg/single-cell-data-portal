import { connect } from "react-redux";
import React from "react";
import * as d3 from "d3";

import { Classes } from "@blueprintjs/core";
import * as globals from "../../../globals";
// @ts-expect-error ts-migrate(2307) FIXME: Cannot find module '../categorical.css' or its cor... Remove this comment to see the full error message
import styles from "../categorical.css";
import Truncate from "../../util/truncate";

import { AnnotationsHelpers } from "../../../util/stateManager";
import { labelPrompt, isLabelErroneous } from "../labelUtil";
import actions from "../../../actions";
import MiniHistogram from "../../miniHistogram";
import MiniStackedBar from "../../miniStackedBar";
import { Dataframe, ContinuousHistogram } from "../../../util/dataframe";
import { track } from "../../../analytics";
import { EVENTS } from "../../../analytics/events";
import { RootState, AppDispatch } from "../../../reducers";
import { Schema, Category } from "../../../common/types/schema";
import { isDataframeDictEncodedColumn } from "../../../util/dataframe/types";
import { AnnotationsState } from "../../../reducers/annotations";
import { CategorySummary } from "../../../util/stateManager/controlsHelpers";
import { ColorTable } from "../../../util/stateManager/colorHelpers";

const STACKED_BAR_HEIGHT = 11;
const STACKED_BAR_WIDTH = 100;

function _currentLabelAsString(label: Category): string {
  return String(label);
}
interface PureCategoryValueProps {
  metadataField: string;
  colorMode: string;
  categoryIndex: number;
  categorySummary: CategorySummary;
  colorAccessor: string;
  colorTable: ColorTable;
  colorData: Dataframe | null;
  categoryData: Dataframe;
  isUserAnno: boolean;
}

interface StateProps {
  annotations: AnnotationsState;
  schema: Schema;
  isDilated: boolean;
  isSelected: boolean;
  label: string;
  labelName: string;
  isColorBy: boolean;
}

interface DispatchProps {
  dispatch: AppDispatch;
}

type Props = StateProps & PureCategoryValueProps & DispatchProps;

const mapDispatchToProps = (dispatch: AppDispatch): DispatchProps => ({
  dispatch,
});

const mapStateToProps = (
  state: RootState,
  ownProps: PureCategoryValueProps,
): StateProps => {
  const { pointDilation, categoricalSelection } = state;
  const {
    metadataField,
    categorySummary,
    colorAccessor,
    colorMode,
    categoryIndex,
    categoryData,
  } = ownProps;

  const label = categorySummary.categoryValues[categoryIndex];
  const isDilated =
    pointDilation.metadataField === metadataField &&
    pointDilation.categoryField === _currentLabelAsString(label);

  const category = categoricalSelection[metadataField];
  const col = categoryData.icol(0);
  const labelName = isDataframeDictEncodedColumn(col)
    ? col.codeMapping[parseInt(label as string, 10)]
    : label;
  const isSelected = category.get(label) ?? true;

  const isColorBy =
    metadataField === colorAccessor &&
    colorMode === "color by categorical metadata";
  return {
    annotations: state.annotations,
    schema: state.annoMatrix?.schema,
    isDilated,
    isSelected,
    label: label as string,
    labelName: labelName as string,
    isColorBy,
  };
};
interface InternalStateProps {
  editedLabelText: string;
}
class CategoryValue extends React.Component<Props, InternalStateProps> {
  constructor(props: Props) {
    super(props);
    this.state = {
      editedLabelText: this.currentLabelAsString(),
    };
  }

  componentDidUpdate(prevProps: Props): void {
    const { metadataField, categoryIndex, categorySummary, colorMode } =
      this.props;
    if (
      prevProps.metadataField !== metadataField ||
      prevProps.colorMode !== colorMode ||
      prevProps.categoryIndex !== categoryIndex ||
      prevProps.categorySummary !== categorySummary
    ) {
      this.setState({
        editedLabelText: this.currentLabelAsString(),
      });
    }
  }

  labelNameError = (name: any) => {
    const { metadataField, schema } = this.props;
    if (name === this.currentLabelAsString()) return false;
    return isLabelErroneous(name, metadataField, schema);
  };

  instruction = (label: any) =>
    labelPrompt(this.labelNameError(label), "New, unique label", ":");

  activateEditLabelMode = () => {
    const { dispatch, metadataField, categoryIndex, label } = this.props;
    dispatch({
      type: "annotation: activate edit label mode",
      metadataField,
      categoryIndex,
      label,
    });
  };

  cancelEditMode = () => {
    const { dispatch, metadataField, categoryIndex, label } = this.props;
    this.setState({
      editedLabelText: this.currentLabelAsString(),
    });
    dispatch({
      type: "annotation: cancel edit label mode",
      metadataField,
      categoryIndex,
      label,
    });
  };

  toggleOff = () => {
    track(EVENTS.EXPLORER_CATEGORICAL_VALUE_SELECT_BUTTON_CLICKED);

    const { dispatch, metadataField, categoryIndex, categorySummary } =
      this.props;

    const label = categorySummary.categoryValues[categoryIndex];
    dispatch(
      actions.selectCategoricalMetadataAction(
        "categorical metadata filter deselect",
        metadataField,
        categorySummary.allCategoryValues,
        label,
        false,
      ),
    );
  };

  shouldComponentUpdate = (
    nextProps: Props,
    nextState: InternalStateProps,
  ): boolean => {
    /*
    Checks to see if at least one of the following changed:
    * world state
    * the color accessor (what is currently being colored by)
    * if this categorical value's selection status has changed
    * the crossfilter (ie, global selection state)
    * the color mode (type of coloring occurring)

    If and only if true, update the component
    */
    const { state } = this;
    const { props } = this;
    const { categoryIndex, categorySummary, isSelected } = props;
    const {
      categoryIndex: newCategoryIndex,
      categorySummary: newCategorySummary,
      isSelected: newIsSelected,
    } = nextProps;

    const label = categorySummary.categoryValues[categoryIndex];
    const newLabel = newCategorySummary.categoryValues[newCategoryIndex];
    const labelChanged = label !== newLabel;
    const valueSelectionChange = isSelected !== newIsSelected;

    const colorAccessorChange = props.colorAccessor !== nextProps.colorAccessor;
    const annotationsChange = props.annotations !== nextProps.annotations;
    const colorModeChange = props.colorMode !== nextProps.colorMode;
    const editingLabel = state.editedLabelText !== nextState.editedLabelText;
    const dilationChange = props.isDilated !== nextProps.isDilated;

    const count = categorySummary.categoryValueCounts[categoryIndex];
    const newCount = newCategorySummary.categoryValueCounts[newCategoryIndex];
    const countChanged = count !== newCount;

    // If the user edits an annotation that is currently colored-by, colors may be re-assigned.
    // This test is conservative - it may cause re-rendering of entire category (all labels)
    // if any one changes, but only for the currently colored-by category.
    const colorMightHaveChanged =
      nextProps.colorAccessor === nextProps.metadataField &&
      props.categorySummary !== nextProps.categorySummary;

    return (
      labelChanged ||
      valueSelectionChange ||
      colorAccessorChange ||
      annotationsChange ||
      editingLabel ||
      dilationChange ||
      countChanged ||
      colorMightHaveChanged ||
      colorModeChange
    );
  };

  toggleOn = () => {
    track(EVENTS.EXPLORER_CATEGORICAL_VALUE_SELECT_BUTTON_CLICKED);

    const { dispatch, metadataField, categoryIndex, categorySummary } =
      this.props;
    const label = categorySummary.categoryValues[categoryIndex];
    dispatch(
      actions.selectCategoricalMetadataAction(
        "categorical metadata filter select",
        metadataField,
        categorySummary.allCategoryValues,
        label,
        true,
      ),
    );
  };

  handleMouseEnter = () => {
    const { dispatch, metadataField, categoryIndex, label } = this.props;
    dispatch({
      type: "category value mouse hover start",
      metadataField,
      categoryIndex,
      label,
    });
  };

  handleMouseExit = () => {
    const { dispatch, metadataField, categoryIndex, label } = this.props;
    dispatch({
      type: "category value mouse hover end",
      metadataField,
      categoryIndex,
      label,
    });
  };

  handleTextChange = (text: string) => {
    this.setState({ editedLabelText: text });
  };

  createHistogramBins = (
    metadataField: string,
    categoryData: Dataframe,
    _colorAccessor: string,
    colorData: Dataframe,
    categoryValue: string,
    width: number,
    height: number,
  ) => {
    /*
      Knowing that colorScale is based off continuous data,
      createHistogramBins fetches the continuous data in relation to the cells relevant to the category value.
      It then separates that data into 50 bins for drawing the mini-histogram
    */
    const groupBy = categoryData.col(metadataField);
    const col = colorData.icol(0);
    const range = col.summarizeContinuous();

    const histogramMap = col.histogramContinuousBy(
      50,
      [range.min, range.max],
      groupBy,
    );

    const bins = histogramMap.has(categoryValue)
      ? (histogramMap.get(categoryValue) as ContinuousHistogram)
      : new Array<number>(50).fill(0);

    const xScale = d3.scaleLinear().domain([0, bins.length]).range([0, width]);

    const largestBin = Math.max(...bins);

    const yScale = d3.scaleLinear().domain([0, largestBin]).range([0, height]);

    return {
      xScale,
      yScale,
      bins,
    };
  };

  createStackedGraphBins = (
    metadataField: string,
    categoryData: Dataframe,
    colorAccessor: string,
    colorData: Dataframe,
    categoryValue: string,
    _colorTable: ColorTable,
    schema: Schema,
    width: number,
  ) => {
    /*
      Knowing that the color scale is based off of categorical data,
      createOccupancyStack obtains a map showing the number if cells per colored value
      Using the colorScale a stack of colored bars is drawn representing the map
     */
    const groupBy = categoryData.col(metadataField);
    const occupancyMap = colorData
      .col(colorAccessor)
      .histogramCategoricalBy(groupBy);

    const occupancy = occupancyMap.get(categoryValue);

    if (occupancy && occupancy.size > 0) {
      // not all categories have occupancy, so occupancy may be undefined.
      const scale = d3
        .scaleLinear()
        /* get all the keys d[1] as an array, then find the sum */
        .domain([0, d3.sum(Array.from(occupancy.values()))])
        .range([0, width]);
      const { categories } = schema.annotations.obsByName[colorAccessor];

      const dfColumn = colorData.col(colorAccessor);
      const categoryValues = dfColumn.summarizeCategorical().categories;

      return {
        domainValues: categoryValues,
        scale,
        domain: categories,
        occupancy,
      };
    }
    return null;
  };

  // If coloring by and this isn't the colorAccessor and it isn't being edited
  shouldRenderStackedBarOrHistogram() {
    const { colorAccessor, isColorBy, annotations } = this.props;
    return !!colorAccessor && !isColorBy && !annotations.isEditingLabelName;
  }

  currentLabelAsString() {
    const { labelName } = this.props;
    return _currentLabelAsString(labelName);
  }

  isAddCurrentSelectionDisabled(crossfilter: any, category: any, value: any) {
    /*
    disable "add current selection to label", if one of the following is true:
    1. no cells are selected
    2. all currently selected cells already have this label, on this category
    */
    const { categoryData } = this.props;

    // 1. no cells selected?
    if (crossfilter.countSelected() === 0) {
      return true;
    }
    // 2. all selected cells already have the label
    const mask = crossfilter.allSelectedMask();
    if (
      AnnotationsHelpers.allHaveLabelByMask(categoryData, category, value, mask)
    ) {
      return true;
    }
    // else, don't disable
    return false;
  }

  renderMiniStackedBar = (): JSX.Element | null => {
    const {
      colorAccessor,
      metadataField,
      categoryData,
      colorData,
      colorTable,
      schema,
      label,
      colorMode,
    } = this.props;
    const isColorBy =
      metadataField === colorAccessor &&
      colorMode === "color by categorical metadata";

    if (!schema) return null;
    if (
      !this.shouldRenderStackedBarOrHistogram ||
      colorMode === "color by expression" ||
      !AnnotationsHelpers.isCategoricalAnnotation(schema, colorAccessor) ||
      isColorBy ||
      !colorData
    ) {
      return null;
    }

    const { domainValues, scale, domain, occupancy } =
      this.createStackedGraphBins(
        metadataField,
        categoryData,
        colorAccessor,
        colorData,
        label,
        colorTable,
        schema,
        STACKED_BAR_WIDTH,
      ) ?? {};

    if (!domainValues || !scale || !domain || !occupancy) {
      return null;
    }

    return (
      <MiniStackedBar
        {...{
          colorTable,
          domainValues,
          scale,
          domain,
          occupancy,
        }}
        // @ts-expect-error ts-migrate(2322) FIXME: Type '{ height: number; width: number; colorTable:... Remove this comment to see the full error message
        height={STACKED_BAR_HEIGHT}
        width={STACKED_BAR_WIDTH}
      />
    );
  };

  renderMiniHistogram = (): JSX.Element | null => {
    const {
      colorAccessor,
      metadataField,
      colorData,
      categoryData,
      colorTable,
      schema,
      label,
      colorMode,
    } = this.props;
    const colorScale = colorTable?.scale;
    if (!schema) return null;
    if (
      !this.shouldRenderStackedBarOrHistogram ||
      // This function returns true on categorical annotations(when stacked bar should not render),
      //  in cases where the colorAccessor is a gene this function will return undefined since genes do not live on the schema
      (AnnotationsHelpers.isCategoricalAnnotation(schema, colorAccessor) ===
        true &&
        colorMode !== "color by expression" &&
        colorMode !== "color by continuous metadata") ||
      !colorData
    ) {
      return null;
    }

    const { xScale, yScale, bins } =
      this.createHistogramBins(
        metadataField,
        categoryData,
        colorAccessor,
        colorData,
        label,
        STACKED_BAR_WIDTH,
        STACKED_BAR_HEIGHT,
      ) ?? {}; // if createHistogramBins returns empty object assign null to deconstructed

    if (!xScale || !yScale || !bins) return null;
    return (
      <MiniHistogram
        {...{
          colorScale,
          xScale,
          yScale,
          bins,
        }}
        // @ts-expect-error ts-migrate(2322) FIXME: Type '{ obsOrVarContinuousFieldDisplayName: any; d... Remove this comment to see the full error message
        obsOrVarContinuousFieldDisplayName={colorAccessor}
        domainLabel={this.currentLabelAsString()}
        height={STACKED_BAR_HEIGHT}
        width={STACKED_BAR_WIDTH}
      />
    );
  };

  render(): JSX.Element {
    const {
      metadataField,
      categoryIndex,
      colorAccessor,
      colorTable,
      isUserAnno,
      isDilated,
      isSelected,
      categorySummary,
      label,
      colorMode,
    } = this.props;
    const colorScale = colorTable?.scale;

    const count = categorySummary.categoryValueCounts[categoryIndex];
    const displayString = this.currentLabelAsString();

    /* this is the color scale, so add swatches below */
    const isColorBy =
      metadataField === colorAccessor &&
      colorMode === "color by categorical metadata";
    const { categoryValueIndices } = categorySummary;

    const valueToggleLabel = `value-toggle-checkbox-${metadataField}-${displayString}`;

    const LEFT_MARGIN = 60;
    const CHECKBOX = 26;
    const CELL_NUMBER = 50;
    const ANNO_MENU = 26;
    const LABEL_MARGIN = 16;
    const CHART_MARGIN = 24;

    const otherElementsWidth =
      LEFT_MARGIN +
      CHECKBOX +
      CELL_NUMBER +
      LABEL_MARGIN +
      (isUserAnno ? ANNO_MENU : 0);

    const labelWidth =
      colorAccessor && !isColorBy
        ? globals.leftSidebarWidth -
          otherElementsWidth -
          STACKED_BAR_WIDTH -
          CHART_MARGIN
        : globals.leftSidebarWidth - otherElementsWidth;

    return (
      <div
        className={
          /* This code is to change the styles on centroid label hover is causing over-rendering */
          `${styles.value}${isDilated ? ` ${styles.hover}` : ""}`
        }
        data-testclass="categorical-row"
        style={{
          padding: "4px 0px 4px 7px",
          display: "flex",
          alignItems: "baseline",
          justifyContent: "space-between",
          marginBottom: "2px",
          borderRadius: "2px",
        }}
        onMouseEnter={this.handleMouseEnter}
        onMouseLeave={this.handleMouseExit}
      >
        <div
          style={{
            margin: 0,
            padding: 0,
            userSelect: "none",
            width: globals.leftSidebarWidth - 145,
            display: "flex",
            justifyContent: "space-between",
          }}
        >
          <div style={{ display: "flex", alignItems: "baseline" }}>
            <label
              htmlFor={valueToggleLabel}
              className={`${Classes.CONTROL} ${Classes.CHECKBOX}`}
              style={{ margin: 0 }}
            >
              <input
                id={valueToggleLabel}
                onChange={isSelected ? this.toggleOff : this.toggleOn}
                data-testclass="categorical-value-select"
                data-testid={`categorical-value-select-${metadataField}-${displayString}`}
                checked={isSelected}
                type="checkbox"
              />
              <span
                className={Classes.CONTROL_INDICATOR}
                onMouseEnter={this.handleMouseExit}
                onMouseLeave={this.handleMouseEnter}
              />
            </label>
            <Truncate>
              <span
                data-testid={`categorical-value-${metadataField}-${displayString}`}
                data-testclass="categorical-value"
                tabIndex={-1}
                style={{
                  width: labelWidth,
                  color:
                    displayString === globals.unassignedCategoryLabel
                      ? "#ababab"
                      : "black",
                  fontStyle:
                    displayString === globals.unassignedCategoryLabel
                      ? "italic"
                      : "normal",
                  display: "inline-block",
                  overflow: "hidden",
                  lineHeight: "1.1em",
                  height: "1.1em",
                  verticalAlign: "middle",
                  marginRight: LABEL_MARGIN,
                }}
              >
                {displayString}
              </span>
            </Truncate>
          </div>
          <span style={{ flexShrink: 0 }}>
            {this.renderMiniStackedBar()}
            {this.renderMiniHistogram()}
          </span>
        </div>
        <div>
          <span>
            <span
              data-testclass="categorical-value-count"
              data-testid={`categorical-value-count-${metadataField}-${displayString}`}
              style={{
                color:
                  displayString === globals.unassignedCategoryLabel
                    ? "#ababab"
                    : "black",
                fontStyle:
                  displayString === globals.unassignedCategoryLabel
                    ? "italic"
                    : "auto",
              }}
            >
              {count}
            </span>

            <svg
              display={isColorBy && categoryValueIndices ? "auto" : "none"}
              style={{
                marginLeft: 5,
                width: 15,
                height: 15,
                backgroundColor:
                  isColorBy && categoryValueIndices && colorScale
                    ? (colorScale(
                        categoryValueIndices.get(label) ?? 0,
                      ) as string)
                    : "inherit",
              }}
            />
          </span>
        </div>
      </div>
    );
  }
}
export default connect(mapStateToProps, mapDispatchToProps)(CategoryValue);
