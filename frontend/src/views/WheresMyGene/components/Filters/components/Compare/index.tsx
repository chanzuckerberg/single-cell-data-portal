import {
  DefaultDropdownMenuOption,
  InputDropdownProps as IInputDropdownProps,
} from "czifui";
import { useContext, useMemo } from "react";
import {
  DispatchContext,
  State,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import { selectCompare } from "src/views/WheresMyGene/common/store/actions";
import { Wrapper, Label, StyledDropdown } from "../common/style";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { LabelWrapper, NewChip } from "./style";
import { COMPARE_OPTIONS } from "src/views/WheresMyGene/common/constants";

const DEFAULT_INPUT_DROPDOWN_PROPS: Partial<IInputDropdownProps> = {
  sdsStyle: "square",
};

interface Props {
  areFiltersDisabled: boolean;
}

// (alec) This is a hack to prevent the group by event from firing when the page loads
let groupByEventTriggeredOnce = false;

export default function Compare({ areFiltersDisabled }: Props): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const { compare } = useContext(StateContext);

  const InputDropdownProps = useMemo(
    () => ({
      ...DEFAULT_INPUT_DROPDOWN_PROPS,
      disabled: areFiltersDisabled,
    }),
    [areFiltersDisabled]
  );

  const optionLabel = useMemo(() => {
    return COMPARE_OPTIONS.find((option) => option.id === compare);
  }, [compare]);

  return (
    <div>
      <LabelWrapper>
        <Label>Group By</Label>
        <NewChip label="NEW" />
      </LabelWrapper>
      <Wrapper>
        <StyledDropdown
          data-testid="compare-dropdown"
          onChange={handleChange}
          label={optionLabel?.name || "None"}
          options={COMPARE_OPTIONS}
          InputDropdownProps={InputDropdownProps}
          value={optionLabel}
        />
      </Wrapper>
    </div>
  );

  function handleChange(value: DefaultDropdownMenuOption | null): void {
    if (!dispatch || !value) return;
    if (!groupByEventTriggeredOnce) {
      track(EVENTS.WMG_OPTION_SELECT_GROUP_BY, {
        group_by_option: value.name,
      });
    } else {
      groupByEventTriggeredOnce = true;
    }

    dispatch(
      selectCompare(
        (value as typeof COMPARE_OPTIONS[number]).id as State["compare"]
      )
    );
  }
}
