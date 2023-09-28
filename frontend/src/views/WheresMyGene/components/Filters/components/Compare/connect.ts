import {
  DefaultDropdownMenuOption,
  InputDropdownProps as IInputDropdownProps,
} from "@czi-sds/components";
import { useCallback, useContext, useMemo } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { COMPARE_OPTIONS } from "src/views/WheresMyGene/common/constants";
import {
  DispatchContext,
  State,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import { selectCompare } from "src/views/WheresMyGene/common/store/actions";
import { Props } from "./types";

export const useConnect = ({ areFiltersDisabled }: Props) => {
  const dispatch = useContext(DispatchContext);
  const { compare } = useContext(StateContext);

  const DEFAULT_INPUT_DROPDOWN_PROPS: Partial<IInputDropdownProps> = useMemo(
    () => ({
      sdsStyle: "square",
    }),
    []
  );

  const InputDropdownProps = useMemo(
    () => ({
      ...DEFAULT_INPUT_DROPDOWN_PROPS,
      disabled: areFiltersDisabled,
    }),
    [areFiltersDisabled, DEFAULT_INPUT_DROPDOWN_PROPS]
  );

  const optionLabel = useMemo(() => {
    return COMPARE_OPTIONS.find((option) => option.id === compare);
  }, [compare]);

  const handleChange = useCallback(
    (value: DefaultDropdownMenuOption | null) => {
      if (!dispatch || !value) return;
      track(EVENTS.WMG_OPTION_SELECT_GROUP_BY, {
        group_by_option: value.name,
      });
      track(EVENTS.WMG_OPTION_SELECT_GROUP_BY, {
        group_by_option: value.name,
      });
      dispatch(
        selectCompare(
          (value as (typeof COMPARE_OPTIONS)[number]).id as State["compare"]
        )
      );
    },
    [dispatch]
  );
  return {
    handleChange,
    optionLabel,
    InputDropdownProps,
  };
};
