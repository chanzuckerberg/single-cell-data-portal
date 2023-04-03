import styled from "@emotion/styled";
import { Button, fontBodyS, fontBodyXs, fontCapsXxs, getColors } from "czifui";
import Modal from "src/components/common/Modal";

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
  width: 600px !important;
  padding: 10px !important;
  h5 {
    color: black;
    height: 28px !important;
    font-size: 24px !important;
    margin: 0px !important;
  }
  padding: 24px !important;
`;

export const StyledTitle = styled.div`
  ${fontCapsXxs}

  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.gray[500]};
    `;
  }}
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
