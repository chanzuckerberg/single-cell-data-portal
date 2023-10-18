import { CommonThemeProps } from "@czi-sds/components";
import styled from "@emotion/styled";
import { textPrimary, textSecondary } from "src/common/theme";

interface StyledTextProps extends CommonThemeProps {
  isInCorpus: boolean;
}
export const StyledText = styled.text<StyledTextProps>`
  fill: ${(props) =>
    props.isInCorpus ? textPrimary(props) : textSecondary(props)};
`;
