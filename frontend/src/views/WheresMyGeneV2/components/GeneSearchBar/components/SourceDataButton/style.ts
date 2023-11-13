import styled from "@emotion/styled";
import {
  ButtonIcon,
  fontBodyXxs,
  fontBodyXxxs,
  getColors,
} from "@czi-sds/components";
import { Badge } from "@mui/base";
import { error400, fontWeightSemibold } from "src/common/theme";

interface ButtonProps {
  refrenceCountSmall?: boolean;
}

export const BadgeCounter = styled(Badge)`
  background-color: ${error400};
  border-radius: 12px;
  width: 28px;
  height: 16px;
  text-align: center;
  position: relative;
  z-index: 22;
  top: 4px;
  left: 14px;

  .MuiBadge-badge {
    position: relative;
    top: -2px;
    margin: 0px;
    ${fontBodyXxxs}
    font-weight: ${fontWeightSemibold};
    color: #fff;
  }
`;

export const StyledButtonDiv = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
  padding: 0 10px;
`;

export const StyledLabel = styled.div`
  ${fontBodyXxs}

  white-space: nowrap;

  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.gray[500]}
    `;
  }}
`;

export const StyledButtonIcon = styled(ButtonIcon)<ButtonProps>`
  width: 30px;
  height: 30px;
  position: absolute;
  cursor: pointer;
  top: 99px;
  z-index: 1;
`;
