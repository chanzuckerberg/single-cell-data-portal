import styled from "styled-components";
import { Image } from "theme-ui";
import { fontStyle } from "../theme";

interface LogoProps {
  size?: number;
}

export const Logo = styled(Image)`
  width: ${({ size = 200 }: LogoProps) => `${size}px`};
  height: 35px;
`;

export const LogoWrapper = styled.span`
  display: flex;
  align-items: center;
`;

interface CorporaWrapperProps {
  small?: boolean;
}

export const CorporaWrapper = styled.span`
  ${fontStyle}
  color: #333333;
  margin-left: 10px;
  padding-left: 10px;
  border-left: 1px solid black;
  font-style: normal;
  font-weight: bold;
  font-size: ${({ small = false }: CorporaWrapperProps) => {
    return small ? "14px" : "24px";
  }};
  line-height: ${({ small = false }: CorporaWrapperProps) => {
    return small ? "16px" : "28px";
  }};
`;
