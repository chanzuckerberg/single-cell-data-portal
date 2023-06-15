import { ECharts, ElementEvent, init } from "echarts";
import {
  memo,
  useEffect,
  useRef,
  useState,
  forwardRef,
  ForwardedRef,
  HTMLAttributes,
} from "react";

import { EMPTY_OBJECT } from "src/common/constants/utils";

import { ChartContainer } from "./style";
import { useUpdateChart } from "./hooks/useUpdateChart";
import { CreateChartOptionsProps } from "./hooks/utils";

export interface Props
  extends HTMLAttributes<HTMLDivElement>,
    CreateChartOptionsProps {
  echartsRendererMode?: "svg" | "canvas";
  onChartMouseMove?: (event: ElementEvent, chart: ECharts) => void;
}

function Chart(props: Props, ref: ForwardedRef<HTMLDivElement>): JSX.Element {
  const {
    width,
    height,
    echartsRendererMode = "canvas",
    onChartMouseMove,
    xAxisData,
    yAxisData,
    data,
    encode,
    itemStyle,
    symbolSize,
    grid,
    ...rest
  } = props;

  if (!width || !height) {
    throw Error("Chart must have width and height > 0");
  }

  const [isChartInitialized, setIsChartInitialized] = useState(false);

  const [chart, setChart] = useState<ECharts | null>(null);
  const innerRef = useRef<HTMLDivElement | null>(null);

  // Initialize charts
  useEffect(() => {
    const { current } = innerRef;

    if (
      !current ||
      isChartInitialized ||
      // (thuang): echart's `init()` will throw error if the container has 0 width or height
      current?.getAttribute("height") === "0" ||
      current?.getAttribute("width") === "0"
    ) {
      return;
    }

    setIsChartInitialized(true);

    const chart = init(current, EMPTY_OBJECT, {
      renderer: echartsRendererMode,
      useDirtyRect: true,
    });

    if (onChartMouseMove) {
      chart.getZr().on("mousemove", function (event: ElementEvent) {
        onChartMouseMove(event, chart);
      });
    }

    setChart(chart);
  }, [innerRef, isChartInitialized, echartsRendererMode, onChartMouseMove]);

  useUpdateChart({
    chart,
    data,
    encode,
    height,
    width,
    xAxisData,
    yAxisData,
    itemStyle,
    symbolSize,
    grid,
  });

  return (
    <ChartContainer height={height} width={width} ref={handleRef} {...rest} />
  );

  function handleRef(element: HTMLDivElement | null) {
    innerRef.current = element;

    if (!ref) return;

    // (thuang): `ref` from `forwardRef` can be a function or a ref object
    if (typeof ref === "function") {
      ref(element);
    } else {
      ref.current = element;
    }
  }
}

export default memo(forwardRef(Chart));
