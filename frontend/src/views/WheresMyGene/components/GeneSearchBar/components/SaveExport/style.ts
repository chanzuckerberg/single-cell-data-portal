import styled from "@emotion/styled";
import FormControlLabel from "@mui/material/FormControlLabel";
import {
  Button,
  fontBodyS,
  fontCapsXxs,
  getColors,
  InputCheckbox,
} from "czifui";
import Modal from "src/components/common/Modal";

export const DownloadButton = styled(Button)`
  ${fontBodyS}
  margin-left: 16px;
  min-width: unset;
  text-transform: unset;
`;
export const StyledDiv = styled.div`
  display: flex;
  flex-direction: column;
  margin: 0;
`;

export const ButtonWrapper = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
  padding: 0 15px;
`;

export const StyledModal = styled(Modal)`
  /* Overwriting some styles for the modal */
  div:first-child {
    padding: 0;
  }
  h5 {
    color: black;
    height: 28px !important;
    font-size: 24px !important;
    margin: 0px !important;
  }
  padding: 24px !important;
  width: 400px !important;
  height: 272px !important;
  min-width: unset !important;
`;

export const StyledTitle = styled.div`
  ${fontCapsXxs}

  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.gray[500]};
    `;
  }}

  margin-bottom: 8px;
`;

export const StyledSection = styled.section`
  margin-top: 16px;
  margin-bottom: 0px;
`;

export const StyledFormControlLabel = styled(FormControlLabel)`
  margin: unset;
  margin-bottom: 8px;
`;

export const StyledInputCheckBox = styled(InputCheckbox)`
  height: 16px;
  width: 16px;
  margin-right: 8px;
`;

export const StyledButtonContainer = styled.div`
  text-align: right;
`;
