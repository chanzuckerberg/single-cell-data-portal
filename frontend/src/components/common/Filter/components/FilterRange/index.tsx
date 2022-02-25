import { Mark, Slider } from "@material-ui/core";
import React, { ChangeEvent, useRef, useState } from "react";
import {
  OnFilterFn,
  Range,
  RangeCategoryView,
} from "src/components/common/Filter/common/entities";
import { formatNumberToMagnitude } from "src/components/common/Filter/common/utils";
import { useSliderStyles } from "src/components/common/Filter/components/FilterRange/style";

type SliderEl = HTMLSpanElement; /* TODO(cc) review */
type SliderRange = number | number[];

interface Props {
  categoryView: RangeCategoryView;
  onFilter: OnFilterFn;
}

/**
 * Returns marks that indicate predetermined values to which the user can move the slider.
 * @param values
 * @returns Array of slider marks.
 */
function buildMarks(values: number[]): Mark[] {
  return values.map((value) => {
    return { label: formatNumberToMagnitude(value), value: value };
  });
}

export default function FilterRange({
  categoryView,
  onFilter,
}: Props): JSX.Element {
  const { key, max, min, selectedMax, selectedMin } = categoryView;
  const classes = useSliderStyles();
  const sliderRef = useRef<SliderEl>(null);
  const [range, setRange] = useState<SliderRange>([
    selectedMin || min,
    selectedMax || max,
  ]);
  const marks = buildMarks([min, max]);

  const onChangeSliderRange = (
    // eslint-disable-next-line @typescript-eslint/ban-types
    _changeEvent: ChangeEvent<{}>, // TODO(cc)
    newRange: SliderRange
  ): void => {
    setRange(newRange);
  };

  const onCommittedSliderRange = (
    // eslint-disable-next-line @typescript-eslint/ban-types
    _changeEvent: ChangeEvent<{}>, // TODO(cc)
    committedRange: SliderRange
  ): void => {
    onFilter(key, committedRange as Range);
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
      ref={sliderRef} // TODO(cc) for wrapper component
      value={range}
      valueLabelFormat={formatNumberToMagnitude}
      valueLabelDisplay="on"
    />
  );
}
