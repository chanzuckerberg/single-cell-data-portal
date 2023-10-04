import {
  DefaultDropdownMenuOption,
  InputDropdownProps as IInputDropdownProps,
} from "@czi-sds/components";
import { useCallback, useContext, useMemo, useRef } from "react";
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

const DEFAULT_INPUT_DROPDOWN_PROPS: Partial<IInputDropdownProps> = {
  sdsStyle: "square",
};

export const useConnect = ({ areFiltersDisabled }: Props) => {
  const dispatch = useContext(DispatchContext);
  const { compare } = useContext(StateContext);

  const InputDropdownProps = useMemo(
    () => ({
      ...DEFAULT_INPUT_DROPDOWN_PROPS,
      disabled: areFiltersDisabled,
    }),
    [areFiltersDisabled]
  );

  const optionLabel: DefaultDropdownMenuOption | undefined = useMemo(() => {
    return COMPARE_OPTIONS.find((option) => option.id === compare);
  }, [compare]);

  const optionLabelRef = useRef<DefaultDropdownMenuOption | undefined>(
    optionLabel
  );

  /**
   * TEMP FIX
   * Warning(thuang): `handleChange` CANNOT depend on `optionLabel`, since a new
   * handleChange passed to SDS Dropdown will trigger another onChange event with the old value,
   * causing infinite loop
   */
  const handleChange = useCallback(
    (value: DefaultDropdownMenuOption | null) => {
      if (!dispatch || !value || optionLabelRef.current === value) return;

      track(EVENTS.WMG_OPTION_SELECT_GROUP_BY, {
        group_by_option: value.name,
      });

      optionLabelRef.current = value;

      dispatch(
        selectCompare(
          (value as (typeof COMPARE_OPTIONS)[number]).id as State["compare"]
        )
      );
    },
    [dispatch, optionLabelRef]
  );
  return {
    handleChange,
    optionLabel,
    InputDropdownProps,
  };
};
