import Slider from "@mui/material/Slider";
import React, { SyntheticEvent, useState } from "react";
import {
  OnFilterFn,
  RangeCategoryView,
} from "src/components/common/Filter/common/entities";
import { formatNumberToScale } from "src/components/common/Filter/common/utils";
import { useSliderStyles } from "src/components/common/Filter/components/FilterContent/components/FilterRange/style";

/**
 * Value returned on change events from MUI Slider.
 */
type SliderRange = number | number[];

/**
 * Slider step size - default to 100.
 */
const STEP_SIZE = 100;

interface Props {
  categoryView: RangeCategoryView;
  onFilter: OnFilterFn;
}

/**
 * Returns marks that indicate predetermined values to which the user can move the slider.
 * @param values
 * @returns Array of slider marks.
 */
function buildMarks(values: number[]) {
  return values.map((value) => {
    return { label: formatNumberToScale(value), value: value };
  });
}

export default function FilterRange({
  categoryView,
  onFilter,
}: Props): JSX.Element {
  const { categoryFilterId, max, min, selectedMax, selectedMin } = categoryView;
  const classes = useSliderStyles();
  const [range, setRange] = useState<SliderRange>([
    selectedMin || min,
    selectedMax || max,
  ]);
  const marks = buildMarks([min, max]);
  const step = Math.floor((max - min) / STEP_SIZE);

  const onChangeSliderRange = (
    _changeEvent: Event,
    newRange: SliderRange
  ): void => {
    setRange(newRange);
  };

  const onCommittedSliderRange = (
    _changeEvent: SyntheticEvent | Event,
    committedRange: SliderRange
  ): void => {
    const [min, max] = committedRange as number[]; // Always expecting a min/max array here.
    onFilter(categoryFilterId, [min, max], [min, max]);
  };

  return (
    <Slider
      className={classes.root}
      color="primary"
      marks={marks}
      max={max}
      min={min}
      onChange={onChangeSliderRange}
      onChangeCommitted={onCommittedSliderRange}
      step={step}
      value={range}
      valueLabelFormat={formatNumberToScale}
      valueLabelDisplay="on"
    />
  );
}
