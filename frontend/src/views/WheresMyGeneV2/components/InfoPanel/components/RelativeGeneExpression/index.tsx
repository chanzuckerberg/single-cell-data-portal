import Image from "next/image";
import { Tooltip } from "@czi-sds/components";
import { Content, Label, LowHigh } from "../../common/style";
import plasmaImage from "./plasma.png";
import { ContentWrapper, Wrapper, TooltipLink } from "./style";
import { MAX_EXPRESSION_LABEL_TEST_ID } from "./constants";
import {
  StyledIconImage,
  StyledTooltip,
  TooltipContent,
  TooltipButton,
} from "src/views/WheresMyGeneV2/components/CellInfoSideBar/style";
import questionMarkIcon from "src/common/images/question-mark-icon.svg";
import { COLOR_LEGEND, COLOR_LEGEND_TOOLTIP_CONTENT } from "./constants";
import { ROUTES } from "src/common/constants/routes";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { useConnect } from "../Legend/connect";

const CONTENT_WIDTH_PX = 120;

interface Props {
  isScaled: boolean;
  maxExpression: number;
}
export default function RelativeGeneExpression({
  isScaled,
  maxExpression,
}: Props): JSX.Element {
  const { setHoverStartTime, handleLegendHoverEnd } = useConnect();
  return (
    <Wrapper id="relative-gene-expression">
      <Label id="relative-gene-expression-label">Gene Expression</Label>
      <Tooltip
        sdsStyle="dark"
        placement="bottom-end"
        width="wide"
        className="fmg-tooltip-icon"
        arrow
        onOpen={() => setHoverStartTime(Date.now())}
        onClose={handleLegendHoverEnd}
        title={
          <StyledTooltip>
            <TooltipContent>
              {COLOR_LEGEND_TOOLTIP_CONTENT}
              <TooltipLink
                href={ROUTES.WMG_DOCS}
                rel="noopener"
                target="_blank"
                onClick={() => {
                  track(EVENTS.WMG_FMG_DOCUMENTATION_CLICKED, {
                    label: COLOR_LEGEND,
                  });
                }}
              >
                our documentation.
              </TooltipLink>
            </TooltipContent>
          </StyledTooltip>
        }
      >
        <TooltipButton sdsStyle="minimal" sdsType="secondary" isAllCaps={false}>
          <StyledIconImage src={questionMarkIcon} />
        </TooltipButton>
      </Tooltip>
      <ContentWrapper>
        <Content>
          <Image
            id="visualization-color-scale"
            src={plasmaImage}
            alt="visualization color scale: interpolateMagma(1.0 - meanExpression)"
            width={CONTENT_WIDTH_PX}
          />
          {(isScaled || maxExpression !== -Infinity) && (
            <LowHigh className="low-high">
              <span>0.0</span>
              <span data-testid={MAX_EXPRESSION_LABEL_TEST_ID}>
                {isScaled ? "1.0" : maxExpression.toFixed(2)}
              </span>
            </LowHigh>
          )}
        </Content>
      </ContentWrapper>
    </Wrapper>
  );
}
