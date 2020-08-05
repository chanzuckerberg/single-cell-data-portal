import styled, { css } from "styled-components";
import { fontStyle as rawFontStyle } from "../common/theme";

const fontStyle = css`
  ${rawFontStyle}
  font-weight: normal;
  line-height: 22px;
`;

export const Wrapper = styled.div`
  ${fontStyle}
  background-color: white;
  position: fixed;
  bottom: 32px;
  left: 32px;
  width: 550px;
  height: 185px;
  border-radius: 4px;
  padding: 32px;
  box-shadow: 0px 5px 9px 1px rgba(169, 169, 169, 0.75);

  letter-spacing: 0.15px;
`;

export const Link = styled.a`
  ${fontStyle}
  text-decoration: none;
`;

export const ButtonWrapper = styled.div`
  margin-top: 16px;
`;

const sharedButtonStyle = css`
  ${rawFontStyle}
  border-radius: 3px;
  font-weight: 600;
  line-height: 16px;
  align-items: center;
  text-align: center;
  padding: 7px;
  border-style: none;
  height: 30px;
  cursor: pointer;
`;

const BUTTON_BLUE = "#0076dc";

export const OKButton = styled.button`
  ${sharedButtonStyle}
  background-color: ${BUTTON_BLUE};
  border-radius: 3px;
  font-weight: 600;
  line-height: 16px;
  align-items: center;
  text-align: center;
  color: #ffffff;
  padding: 7px;
  border-style: none;
  width: 150px;
  height: 30px;
  cursor: pointer;
  margin-right: 10px;
`;

export const NoButton = styled.button`
  ${sharedButtonStyle}
  width: 90px;
  background-color: white;
  color: ${BUTTON_BLUE};
  border: 1px solid ${BUTTON_BLUE};
`;
