import { Tooltip } from "@czi-sds/components";
import questionMarkIcon from "src/common/images/question-mark-icon.svg";
import { StyledDropdown, Wrapper } from "../common/style";
import { Label } from "src/views/WheresMyGeneV2/components/InfoPanel/common/style";
import { LabelWrapper } from "./style";
import {
  StyledIconImage,
  TooltipButton,
  StyledTooltip,
} from "src/views/WheresMyGeneV2/components/CellInfoSideBar/style";
import { COLOR_SCALE_TOOLTIP_TEXT } from "src/views/WheresMyGeneV2/common/constants";
import { COLOR_SCALE_OPTIONS } from "./constants";
import { DEFAULT_INPUT_DROPDOWN_PROPS, Props } from "./types";
import { useConnect } from "./connect";

export default function ColorScale({ setIsScaled }: Props): JSX.Element {
  const { colorScaleOnChange, colorScaledOption } = useConnect({ setIsScaled });

  return (
    <Wrapper>
      <LabelWrapper>
        <Label>Color Scale</Label>
        <Tooltip
          sdsStyle="dark"
          placement="right"
          width="default"
          arrow
          title={
            <StyledTooltip>
              <div>{COLOR_SCALE_TOOLTIP_TEXT}</div>
            </StyledTooltip>
          }
        >
          <TooltipButton
            data-testid="color-scale-tooltip-icon"
            sdsStyle="minimal"
            sdsType="secondary"
            isAllCaps={false}
          >
            <StyledIconImage alt="question mark" src={questionMarkIcon} />
          </TooltipButton>
        </Tooltip>
      </LabelWrapper>

      <StyledDropdown
        data-testid="color-scale-dropdown"
        onChange={colorScaleOnChange}
        label={colorScaledOption.name}
        options={COLOR_SCALE_OPTIONS}
        InputDropdownProps={DEFAULT_INPUT_DROPDOWN_PROPS}
        value={colorScaledOption}
      />
    </Wrapper>
  );
}
