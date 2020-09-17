import { BLUE, GREY } from "src/components/common/theme";
import styled, { css } from "styled-components";

export const Wrapper = styled.div`
  width: 625px;
  height: 300px;
`;

const buttonStyle = css`
  margin-right: 10px;
  border: none;
  cursor: pointer;
`;

export const Cancel = styled.button`
  ${buttonStyle}
  background-color: transparent;
  color: ${GREY.DARK};
  font-size: 14px;
`;

const sharedDownloadStyle = css`
  ${buttonStyle}
  background-color: ${BLUE};
  color: white;
  font-size: 14px;
  width: 80px;
  height: 37px;
  border-radius: 4px;
  margin-right: 10px;
`;

export const DisabledDownload = styled.button`
  ${sharedDownloadStyle}

  cursor: default;
  opacity: 0.6;
  filter: unset;
`;

export const Download = styled.a`
  ${sharedDownloadStyle}

  display: flex;
  align-items: center;
  justify-content: center;

  :hover {
    filter: brightness(1.1);
    color: white;
    cursor: pointer;
    text-decoration: none;
  }
`;
