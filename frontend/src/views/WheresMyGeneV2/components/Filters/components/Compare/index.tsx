import { Wrapper, Label, StyledDropdown } from "../common/style";
import { CompareWrapper, LabelWrapper } from "./style";
import {
  COMPARE_OPTIONS,
  GROUP_BY_TOOLTIP_TEXT,
  SELECT_TISSUE_GENE_TEXT,
} from "src/views/WheresMyGeneV2/common/constants";
import { DefaultAutocompleteOption, Tooltip } from "@czi-sds/components";
import {
  StyledTooltip,
  TooltipButton,
} from "src/views/WheresMyGeneV2/components/CellInfoSideBar/style";
import { Props } from "./types";
import { useConnect } from "./connect";
import { StyledQuestionMarkIcon } from "src/common/style";

export default function Compare({ areFiltersDisabled }: Props): JSX.Element {
  const { optionLabel, handleChange, InputDropdownProps } = useConnect({
    areFiltersDisabled,
  });

  return (
    <CompareWrapper>
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
              <StyledQuestionMarkIcon />
            </TooltipButton>
          </Tooltip>
        </Label>
      </LabelWrapper>
      <Wrapper>
        <StyledDropdown<DefaultAutocompleteOption, false, false, false>
          data-testid="compare-dropdown"
          onChange={handleChange}
          label={optionLabel?.name || "None"}
          options={COMPARE_OPTIONS}
          InputDropdownProps={InputDropdownProps}
          value={optionLabel}
        />
      </Wrapper>
    </CompareWrapper>
  );
}
