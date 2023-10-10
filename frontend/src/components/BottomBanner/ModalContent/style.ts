import {
  Button,
  InputText,
  fontBodyS,
  fontBodyXxxs,
  fontHeaderXl,
} from "@czi-sds/components";
import styled from "@emotion/styled";
import { error400, gray500 } from "src/common/theme";

export const StyledTitle = styled.div`
  ${fontHeaderXl}

  letter-spacing: -0.019em;
  font-size: 24px !important;
  margin: 0;
  height: auto !important;

  padding-top: 16px;
  padding-bottom: 8px;
`;

export const StyledDescription = styled.div`
  ${fontBodyS}

  letter-spacing: -0.006em;
  padding-bottom: 16px;
`;

export const StyledForm = styled.form`
  display: flex;
  flex-direction: row;
  justify-content: space-between;
  height: 34px;
  margin-bottom: 0;
  align-items: center;
  width: 100%;
`;

export const StyledInputText = styled(InputText)`
  .MuiOutlinedInput-root.Mui-focused .MuiOutlinedInput-notchedOutline {
    border-color: #5826c1 !important;
  }

  flex: 1;
  margin-right: 4px;
  margin-bottom: 0px;
  display: inline-flex;
`;

export const StyledSubmitButton = styled(Button)`
  padding: 6px 12px;
  width: 91px;
  height: 34px;
  background: #8f5aff;
  font-weight: 500;

  :hover {
    background: #5826c1;
  }
`;

export const StyledDisclaimer = styled.div`
  ${fontBodyXxxs}

  letter-spacing: -0.005em;

  /*
   * beta intent does not exist for SDS banner, but the colors do
   * targeting specific id to overwrite style
   */
  color: ${gray500};
`;

export const StyledErrorMessage = styled.div`
  ${fontBodyXxxs}

  letter-spacing: -0.005em;

  align-self: flex-start;

  height: 16px;
  margin-top: 4px;
  margin-bottom: 4px;

  /*
   * beta intent does not exist for SDS banner, but the colors do
   * targeting specific id to overwrite style
   */
  color: ${error400};
`;
