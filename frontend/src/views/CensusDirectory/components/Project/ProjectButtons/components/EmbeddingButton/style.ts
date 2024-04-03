import {
  DialogContent,
  CommonThemeProps,
  fontCapsXxxs,
  fontCapsXxs,
  Callout,
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
  max-height: 250px;
  overflow-y: scroll;

  pre {
    margin: 0;
  }
  button {
    position: absolute;
    padding: 0;
    color: ${primary300};
    right: 8px;
    font-weight: 500;
  }
  & > button {
    top: 8px;
  }

  /* URI Container */
  & > div {
    position: absolute;
    top: ${({ uriTopPosition }) => uriTopPosition}px;
    left: 0;
    width: 100%;
    height: ${({ lineHeight }) => lineHeight}px;
    background: rgba(255, 255, 255, 0.05);
    & > button {
      top: 5px;
    }
  }

  button:hover {
    color: ${primary200};
  }
`;
export const StyledDialogContent = styled(DialogContent)`
  padding: 4px;
  & > * {
    margin: 24px 0;
  }
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
    ${fontCapsXxs}
    font-weight: ${fontWeightSemibold};
    background: #fff;
    content: "OR";
    padding: 0 4px;
    position: relative;
    top: -10px;
  }
`;
export const StyledCallout = styled(Callout)`
  width: 100%;
`;
