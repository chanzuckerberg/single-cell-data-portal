import { Wrapper, Label, StyledDropdown } from "../common/style";
import { LabelWrapper } from "./style";
import {
  COMPARE_OPTIONS,
  GROUP_BY_TOOLTIP_TEXT,
  SELECT_TISSUE_GENE_TEXT,
} from "src/views/WheresMyGeneV2/common/constants";
import { Tooltip } from "@czi-sds/components";
import {
  StyledIconImage,
  StyledTooltip,
  TooltipButton,
} from "src/views/WheresMyGeneV2/components/CellInfoSideBar/style";
import questionMarkIcon from "src/common/images/question-mark-icon.svg";
import { Props } from "./types";
import { useConnect } from "./connect";

export default function Compare({ areFiltersDisabled }: Props): JSX.Element {
  const { optionLabel, handleChange, InputDropdownProps } = useConnect({
    areFiltersDisabled,
  });

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
        {/* <NewChip label="NEW" /> */}
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
}
