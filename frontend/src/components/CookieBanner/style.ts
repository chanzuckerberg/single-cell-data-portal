import styled, { css } from "styled-components";
import { fontStyle as rawFontStyle, OLD_BLUE } from "../common/theme";

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
  text-align: center;
  padding: 7px;
  border-style: none;
  height: 30px;
  cursor: pointer;
`;

export const OKButton = styled.button`
  ${sharedButtonStyle}
  background-color: ${OLD_BLUE};
  color: white;
  margin-right: 10px;
`;

export const NoButton = styled.button`
  ${sharedButtonStyle}
  background-color: white;
  border: 1px solid ${OLD_BLUE};
  color: ${OLD_BLUE};
  width: 90px;
`;
