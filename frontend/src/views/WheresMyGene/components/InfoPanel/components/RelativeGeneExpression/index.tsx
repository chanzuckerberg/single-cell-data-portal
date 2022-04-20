import { Checkbox, Icon, Tooltip } from "czifui";
import Image from "next/image";
import { Content, Header, LowHigh } from "../../common/style";
import plasmaImage from "./plasma.png";
import {
  ContentWrapper,
  LabelWrapper,
  StyledFormControlLabel,
  Wrapper,
} from "./style";

const CONTENT_WIDTH_PX = 100;

interface Props {
  handleIsScaledChange: () => void;
  isScaled: boolean;
}
export default function RelativeGeneExpression({
  handleIsScaledChange,
  isScaled,
}: Props): JSX.Element {
  return (
    <Wrapper>
      <Header>Gene Expression</Header>
      <ContentWrapper>
        <Content>
          <Image
            src={plasmaImage}
            alt="visualization color scale: interpolateMagma(1.0 - meanExpression)"
            width={CONTENT_WIDTH_PX}
          />
          <LowHigh>
            <span>0.0</span>
            <span>{isScaled ? "1.0" : "6.0"}</span>
          </LowHigh>
        </Content>
        <StyledFormControlLabel
          control={
            <Checkbox
              inputProps={{
                "aria-label":
                  "is scaled? " + isScaled ? "checked" : "unchecked",
                checked: isScaled,
                hidden: true,
              }}
              onChange={handleIsScaledChange}
              stage={isScaled ? "checked" : "unchecked"}
            />
          }
          label={
            <LabelWrapper>
              <span>Scaled</span>
              <Tooltip title="Expression is scaled to the range [0,1]. Scaling is done by assigning the minimum value in the current view to 0 and the max is assigned to 1.">
                <span>
                  <Icon sdsIcon="infoCircle" sdsSize="xs" sdsType="static" />
                </span>
              </Tooltip>
            </LabelWrapper>
          }
        />
      </ContentWrapper>
    </Wrapper>
  );
}
