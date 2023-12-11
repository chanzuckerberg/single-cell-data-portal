import {
  DialogContent,
  CommonThemeProps,
  fontCapsXxxs,
  fontCapsXxs,
} from "@czi-sds/components";
import styled from "@emotion/styled";
import {
  fontWeightSemibold,
  gray300,
  gray400,
  primary200,
  primary300,
} from "src/common/theme";

interface CodeSnippetProps extends CommonThemeProps {
  uriTopPosition: number;
  lineHeight: number;
}
export const CodeSnippet = styled.div<CodeSnippetProps>`
  position: relative;
  button {
    position: absolute;
    padding: 0;
    color: ${primary300};
    &:hover {
      color: ${primary200};
    }
    right: 8px;
  }

  & > button {
    top: 8px;
  }

  // URI Container
  & > div {
    position: absolute;
    top: ${({ uriTopPosition }) => uriTopPosition}px;
    left: 0;
    width: 100%;
    height: ${({ lineHeight }) => lineHeight}px;
    background: rgba(255, 255, 255, 0.05);
    & > button {
      top: 5px;
      span {
        content: "Copy URI";
      }
    }
  }
`;
export const StyledDialogContent = styled(DialogContent)`
  padding: 4px;
`;
export const Label = styled.label`
  ${fontCapsXxxs}
  font-weight: ${fontWeightSemibold};
  font-feature-settings:
    "clig" off,
    "liga" off;
  color: ${gray400};
`;
export const Break = styled.hr`
  border: none;
  border-top: 1px solid ${gray300};
  color: ${gray400};
  width: 64px;
  margin: 0 auto;
  overflow: visible;
  text-align: center;
  height: 5px;
  background: #fff;

  &:after {
    ${fontCapsXxs};
    font-weight: ${fontWeightSemibold};
    background: #fff;
    content: "OR";
    padding: 0 4px;
    position: relative;
    top: -10px;
  }
`;
