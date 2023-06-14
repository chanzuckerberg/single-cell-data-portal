import { interpolateMagma } from "d3-scale-chromatic";

import Chart, {
  Props,
} from "src/views/WheresMyGene/components/HeatMap/components/Chart/components/Chart";

import { symbolSize } from "src/views/WheresMyGene/components/HeatMap/components/Chart/utils";

import { Wrapper } from "./style";

const MATRIX_SIZE = 5;

const data: { percentage: number; size: number; x: number; y: number }[] = [];

for (const i of Array(MATRIX_SIZE).keys()) {
  for (const j of Array(MATRIX_SIZE).keys()) {
    data.push({
      percentage: (i * j) / (MATRIX_SIZE - 1) ** 2,
      size: i * j,
      x: i,
      y: j,
    });
  }
}

const axisData = Array.from(Array(MATRIX_SIZE)).map((_, i) => i);

const DEFAULT_PROPS: Props = {
  data,
  xAxisData: axisData,
  yAxisData: axisData,
  encode: { x: "x", y: "y" },
  height: 500,
  width: 500,
  itemStyle: {
    color(props) {
      const { percentage } = props.data as {
        percentage: number;
      };

      return interpolateMagma(percentage);
    },
  },
  symbolSize,
  echartsRendererMode: "svg",
  id: "chart",
};

export default function Test() {
  return (
    <Wrapper data-testid="test">
      <div>
        <h3>Default</h3>
        <Chart {...DEFAULT_PROPS} data-testid="chart" />
      </div>

      <div>
        <h3>Canvas mode</h3>
        <Chart
          {...DEFAULT_PROPS}
          data-testid="canvas-chart"
          echartsRendererMode="canvas"
        />
      </div>
    </Wrapper>
  );
}
