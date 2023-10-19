import {
  InputDropdownProps as IInputDropdownProps,
  Tooltip,
} from "@czi-sds/components";
import { useContext, useMemo } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import questionMarkIcon from "src/common/images/question-mark-icon.svg";
import { selectSortBy } from "src/views/WheresMyGene/common/store/actions";
import { SORT_BY } from "src/views/WheresMyGene/common/types";
import { StyledDropdown, Wrapper } from "../common/style";
import { Label } from "../../../InfoPanel/common/style";
import { LabelWrapper } from "./style";
import {
  StyledIconImage,
  TooltipButton,
  StyledTooltip,
} from "src/views/WheresMyGeneV2/components/CellInfoSideBar/style";
import { COLOR_SCALE_TOOLTIP_TEXT } from "src/views/WheresMyGene/common/constants";

interface Props {
  setIsScaled: (prevIsScaled: boolean) => void;
}

// (ashin-czi): Used by SaveExport SVG generation to recreate the color scale plasma image in an SVG
// Not the best solution but importing SVG files creates conflicts with the Loader
// and can't get the actual content of the SVG file needed for SVG creation
export const PLASMA_SVG_STRING = `
  <rect width="120" height="16" fill="url(#paint0_linear_5151_537460)" />
  <defs>
      <linearGradient id="paint0_linear_5151_537460" x1="120" y1="8" x2="8.0516e-07" y2="8.00001"
          gradientUnits="userSpaceOnUse">
          <stop stop-color="#090720" />
          <stop offset="0.145123" stop-color="#36106B" />
          <stop offset="0.28796" stop-color="#6B1D81" />
          <stop offset="0.436206" stop-color="#9C2E7F" />
          <stop offset="0.582062" stop-color="#D3436E" />
          <stop offset="0.719889" stop-color="#F66E5C" />
          <stop offset="0.864407" stop-color="#FEA973" />
          <stop offset="1" stop-color="#FDE2A3" />
      </linearGradient>
  </defs>
`;

const COLOR_SCALE_OPTIONS = [
  { id: SORT_BY.COLOR_SCALED, name: "Scaled" },
  { id: SORT_BY.COLOR_UNSCALED, name: "Unscaled" },
];

const DEFAULT_INPUT_DROPDOWN_PROPS: Partial<IInputDropdownProps> = {
  sdsStyle: "square",
};

export default function ColorScale({ setIsScaled }: Props): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const { sortBy } = useContext(StateContext);

  const colorScaledOption = useMemo(() => {
    return (
      COLOR_SCALE_OPTIONS.find((option) => option.id === sortBy.scaled) ||
      COLOR_SCALE_OPTIONS[0]
    );
  }, [sortBy]);

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
            <StyledIconImage src={questionMarkIcon} />
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
}
