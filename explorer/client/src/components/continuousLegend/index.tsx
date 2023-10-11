import React from "react";
import { Action } from "redux";
import { connect, shallowEqual } from "react-redux";
import * as d3 from "d3";
import { interpolateCool } from "d3-scale-chromatic";
import Async, { AsyncProps } from "react-async";
import AnnoMatrix from "../../annoMatrix/annoMatrix";
import { AppDispatch, RootState } from "../../reducers";

import {
  createColorTable,
  createColorQuery,
  ColorTable,
  ColorRange,
} from "../../util/stateManager/colorHelpers";
import { ColorsState } from "../../reducers/colors";
import { Genesets } from "../../reducers/genesets";
import * as globals from "../../globals";
// create continuous color legend
const continuous = (selectorId: any, colorScale: any, colorAccessor: any) => {
  const legendHeight = 200;
  const legendWidth = 80;
  const margin = { top: 10, right: 60, bottom: 10, left: 2 };

  const canvas = d3
    .select(selectorId)
    .style("height", `${legendHeight}px`)
    .style("width", `${legendWidth}px`)
    .append("canvas")
    .attr("height", legendHeight - margin.top - margin.bottom)
    .attr("width", 1)
    .style("height", `${legendHeight - margin.top - margin.bottom}px`)
    .style("width", `${legendWidth - margin.left - margin.right}px`)
    .style("position", "absolute")
    .style("top", `${margin.top + 1}px`)
    .style("left", `${margin.left + 1}px`)
    .style(
      "transform",
      "scale(1,-1)",
    ) /* flip it! dark is high value light is low.
    we flip the color scale as well [1, 0] instead of [0, 1] */
    .node();

  const ctx = canvas?.getContext("2d");

  const legendScale = d3
    .scaleLinear()
    .range([1, legendHeight - margin.top - margin.bottom])
    .domain([
      colorScale.domain()[1],
      colorScale.domain()[0],
    ]); /* we flip this to make viridis colors dark if high in the color scale */

  // image data hackery based on http://bl.ocks.org/mbostock/048d21cf747371b11884f75ad896e5a5
  const image = ctx?.createImageData(1, legendHeight);
  if (!image) return;
  d3.range(legendHeight).forEach((i) => {
    const c = d3.rgb(colorScale(legendScale.invert(i)));
    image.data[4 * i] = c.r;
    image.data[4 * i + 1] = c.g;
    image.data[4 * i + 2] = c.b;
    image.data[4 * i + 3] = 255;
  });

  ctx?.putImageData(image, 0, 0);

  // A simpler way to do the above, but possibly slower. keep in mind the legend
  // width is stretched because the width attr of the canvas is 1
  // See http://stackoverflow.com/questions/4899799/whats-the-best-way-to-set-a-single-pixel-in-an-html5-canvas
  /*
  d3.range(legendheight).forEach(function(i) {
    ctx.fillStyle = colorscale(legendscale.invert(i));
    ctx.fillRect(0,i,1,1);
  });
  */

  const legendAxis = d3
    .axisRight(legendScale)
    .ticks(6)
    .tickFormat(
      d3.format(
        legendScale.domain().some((n) => Math.abs(n) >= 10000) ? ".0e" : ",",
      ),
    );

  const svg = d3
    .select(selectorId)
    .append("svg")
    .attr("height", `${legendHeight}px`)
    .attr("width", `${legendWidth}px`)
    .style("position", "absolute")
    .style("left", "0px")
    .style("top", "0px");

  svg
    .append("g")
    .attr("class", "axis")
    .attr(
      "transform",
      `translate(${legendWidth - margin.left - margin.right + 3},${
        margin.top
      })`,
    )
    .call(legendAxis);

  // text label for the y axis
  svg
    .append("text")
    .attr("transform", "rotate(-90)")
    .attr("y", 2)
    .attr("x", 0 - legendHeight / 2)
    .attr("dy", "1em")
    .attr("data-testid", "continuous_legend_color_by_label")
    .attr("aria-label", colorAccessor)
    .style("text-anchor", "middle")
    .style("fill", "white")
    .text(colorAccessor);
};

interface FetchedAsyncProps {
  colorAccessor: string;
  colorScale: ColorTable["scale"];
  range?: ColorRange;
  domainMin: number;
  domainMax: number;
  colorMode: Action["type"];
}

interface StateProps {
  annoMatrix: AnnoMatrix;
  colors: ColorsState;
  genesets: Genesets;
}
interface DispatchProps {
  handleColorSuccess: () => void;
}
const mapStateToProps = (state: RootState): StateProps => ({
  annoMatrix: state.annoMatrix,
  colors: state.colors,
  genesets: state.genesets.genesets,
});
const mapDispatchToProps = (dispatch: AppDispatch): DispatchProps => ({
  handleColorSuccess: () =>
    dispatch({ type: "color by geneset mean expression success" }),
});

type Props = StateProps & DispatchProps;

class ContinuousLegend extends React.Component<Props> {
  static watchAsync(props: any, prevProps: any) {
    return !shallowEqual(props.watchProps, prevProps.watchProps);
  }

  cachedWatchProps: StateProps | null;

  cachedAsyncProps: FetchedAsyncProps | null;

  constructor(props: Props) {
    super(props);
    this.cachedWatchProps = null;
    this.cachedAsyncProps = null;
  }

  fetchAsyncProps = async (
    props: AsyncProps<FetchedAsyncProps | null>,
  ): Promise<FetchedAsyncProps | null> => {
    const { annoMatrix, colors, genesets } = props.watchProps as StateProps;

    if (!colors || !annoMatrix || !colors.colorMode || !colors.colorAccessor)
      return null;
    const { schema } = annoMatrix;
    const { colorMode, colorAccessor, userColors } = colors;
    const colorQuery = createColorQuery(
      colorMode,
      colorAccessor,
      schema,
      genesets,
    );
    const colorDf = colorQuery
      ? await annoMatrix.fetch(...colorQuery, globals.numBinsObsX)
      : null;
    const colorTable = createColorTable(
      colorMode,
      colorAccessor,
      colorDf,
      schema,
      userColors,
    );
    const colorScale = colorTable.scale;
    const range = (colorScale?.range ?? (() => [0, 0])) as ColorRange;
    const [domainMin, domainMax] = colorScale?.domain?.() ?? [0, 0];

    const result: FetchedAsyncProps = {
      colorAccessor,
      colorScale,
      range,
      domainMin,
      domainMax,
      colorMode,
    };
    return result;
  };

  updateContinuousLegend = (asyncProps: FetchedAsyncProps) => {
    const { handleColorSuccess } = this.props;
    const {
      colorAccessor,
      colorScale,
      range,
      domainMin,
      domainMax,
      colorMode,
    } = asyncProps;
    if (
      colorAccessor &&
      colorScale &&
      range &&
      (domainMin ?? 0) < (domainMax ?? 0)
    ) {
      /* fragile! continuous range is 0 to 1, not [#fa4b2c, ...], make this a flag? */
      const r = range();
      if (!(typeof r === "function" && r()[0][0] === "#")) {
        continuous(
          "#continuous_legend",
          d3.scaleSequential(interpolateCool).domain(colorScale.domain()),
          colorAccessor,
        );
      }
    }
    if (colorScale && colorMode === "color by geneset mean expression") {
      handleColorSuccess();
    }
  };

  render() {
    const { annoMatrix, colors, genesets } = this.props;
    return (
      <div
        id="continuous_legend"
        style={{
          position: "absolute",
          left: 8,
          top: 35,
          zIndex: 1,
          pointerEvents: "none",
        }}
      >
        <Async
          watchFn={ContinuousLegend.watchAsync}
          promiseFn={this.fetchAsyncProps}
          watchProps={{
            annoMatrix,
            colors,
            genesets,
          }}
        >
          <Async.Fulfilled>
            {(asyncProps: FetchedAsyncProps) => {
              if (
                !shallowEqual(asyncProps, this.cachedAsyncProps) &&
                asyncProps &&
                asyncProps.colorMode &&
                asyncProps.colorMode !== "color by categorical metadata"
              ) {
                d3.select("#continuous_legend").selectAll("*").remove();
                this.updateContinuousLegend(asyncProps);
              }
              this.cachedAsyncProps = asyncProps;
              // if asyncProps is null (no colorAccessor), remove legend
              // or if colorMode is categorical, remove legend
              if (
                !asyncProps ||
                asyncProps.colorMode === "color by categorical metadata"
              )
                d3.select("#continuous_legend").selectAll("*").remove();
              return null;
            }}
          </Async.Fulfilled>
        </Async>
      </div>
    );
  }
}

export default connect(mapStateToProps, mapDispatchToProps)(ContinuousLegend);
