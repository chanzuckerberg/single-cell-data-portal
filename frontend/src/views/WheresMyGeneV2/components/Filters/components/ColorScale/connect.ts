import { useContext, useMemo } from "react";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGeneV2/common/store";
import { SORT_BY } from "src/views/WheresMyGeneV2/common/types";
import { COLOR_SCALE_OPTIONS } from "./constants";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { selectSortBy } from "src/views/WheresMyGeneV2/common/store/actions";
import { Props } from "./types";

export const useConnect = ({
  setIsScaled,
}: {
  setIsScaled: Props["setIsScaled"];
}) => {
  const dispatch = useContext(DispatchContext);
  const { sortBy } = useContext(StateContext);

  const colorScaledOption = useMemo(() => {
    return (
      COLOR_SCALE_OPTIONS.find((option) => option.id === sortBy.scaled) ||
      COLOR_SCALE_OPTIONS[0]
    );
  }, [sortBy]);

  function colorScaleOnChange(
    value: { id?: SORT_BY; name: string } | null
  ): void {
    if (!dispatch || !value || colorScaledOption.name === value.name) return;
    track(EVENTS.WMG_OPTION_SELECT_COLOR_SCALE, {
      color_scale_view_option: value.name,
    });
    setIsScaled(value.id == SORT_BY.COLOR_SCALED ? true : false);
    dispatch(selectSortBy({ scaled: value.id as SORT_BY }));
  }

  return {
    colorScaleOnChange,
    colorScaledOption,
  };
};
