import styled from "@emotion/styled";
import {
  Button,
  fontBodyS,
  fontBodyXs,
  fontCapsXxs,
  getColors,
} from "@czi-sds/components";
import Modal from "src/components/common/Modal";
import { gray500 } from "src/common/theme";

export const DOWNLOAD_MODAL_WIDTH_PX = 600;
export const DOWNLOAD_MODAL_PADDING = 24;

export const DownloadButton = styled(Button)`
  ${fontBodyS}
  margin-left: 16px;
  min-width: unset;
  text-transform: unset;
`;

export const StyledModal = styled(Modal)`
  /* Overriding some styles for the modal to match figma */
  div:first-child {
    padding: 0;
  }
  min-width: ${DOWNLOAD_MODAL_WIDTH_PX}px !important;
  max-width: ${DOWNLOAD_MODAL_WIDTH_PX}px !important;
  padding: ${DOWNLOAD_MODAL_PADDING}px !important;
  h5 {
    color: black;
    height: 28px !important;
    font-size: 24px !important;
    margin: 0px !important;
  }
`;

export const StyledTitle = styled.div`
  ${fontCapsXxs}

  color: ${gray500};
`;

export const StyledModalContent = styled.div`
  padding-top: 8px;
  display: flex;
  flex-direction: column;
`;

export const StyledSection = styled.section`
  padding: 8px 0 8px 0;
  width: 100%;
  display: flex;
  flex-direction: column;
`;

export const StyledInputCheckboxWrapper = styled.div`
  padding-bottom: 30px;
  width: 100%;
  label {
    width: 100%;
    margin-left: -9px;
    &::after {
      font-size: 13px !important;
    }
  }
  span {
    font-size: unset !important;
    font-size: 14px !important;
  }
`;

export const StyledMessage = styled.div`
  ${fontBodyXs}

  display: flex;
  justify-content: center;

  margin-top: 12px;

  border-radius: 4px;

  ${(props) => {
    const colors = getColors(props);

    return `
      background: ${colors?.gray[100]};
      color: ${colors?.gray[500]};
    `;
  }}
`;

export const ButtonWrapper = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
  padding: 0 15px;
`;

export const StyledButtonContainer = styled.div`
  padding-top: 16px;
  text-align: right;
`;
