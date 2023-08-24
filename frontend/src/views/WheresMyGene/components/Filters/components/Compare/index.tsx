import {
  DefaultDropdownMenuOption,
  InputDropdownProps as IInputDropdownProps,
} from "@czi-sds/components";
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
import {
  COMPARE_OPTIONS,
  GROUP_BY_TOOLTIP_TEXT,
  SELECT_TISSUE_GENE_TEXT,
} from "src/views/WheresMyGene/common/constants";
import { Tooltip } from "@czi-sds/components";
import {
  StyledIconImage,
  StyledTooltip,
  TooltipButton,
} from "../../../CellInfoSideBar/style";
import questionMarkIcon from "src/common/images/question-mark-icon.svg";

const DEFAULT_INPUT_DROPDOWN_PROPS: Partial<IInputDropdownProps> = {
  sdsStyle: "square",
};

interface Props {
  areFiltersDisabled: boolean;
}

// (alec) This is a hack to prevent the group by event from firing when the page loads
let groupByEventTriggeredAtLeastOnce = false;

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
        <Label>
          Group By
          <Tooltip
            className="group-by-tooltip-icon"
            sdsStyle="dark"
            placement="right"
            width="default"
            arrow
            title={
              <StyledTooltip>
                {areFiltersDisabled && (
                  <p>
                    <em>{SELECT_TISSUE_GENE_TEXT}</em>
                  </p>
                )}
                <div>{GROUP_BY_TOOLTIP_TEXT}</div>
              </StyledTooltip>
            }
          >
            <TooltipButton
              data-testid="group-by-tooltip-icon"
              sdsStyle="minimal"
              sdsType="secondary"
              isAllCaps={false}
            >
              <StyledIconImage src={questionMarkIcon} />
            </TooltipButton>
          </Tooltip>
        </Label>
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
    if (groupByEventTriggeredAtLeastOnce) {
      track(EVENTS.WMG_OPTION_SELECT_GROUP_BY, {
        group_by_option: value.name,
      });
    } else if (value.name !== "None") {
      groupByEventTriggeredAtLeastOnce = true;
      track(EVENTS.WMG_OPTION_SELECT_GROUP_BY, {
        group_by_option: value.name,
      });
    }

    dispatch(
      selectCompare(
        (value as (typeof COMPARE_OPTIONS)[number]).id as State["compare"]
      )
    );
  }
}
