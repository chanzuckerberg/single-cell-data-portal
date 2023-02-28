import styled from "@emotion/styled";
import { ButtonIcon, fontBodyXxs, getColors, Notification } from "czifui";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";
import { LegendWrapper } from "../../../InfoPanel/components/Legend/style";

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

export const StyledButtonIcon = styled(ButtonIcon)`
  width: 30px;
  height: 30px;
`;

export const StyledNotification = styled(Notification)`
  ${LegendWrapper} .MuiAlert-root {
    position: absolute;
    top: ${HEADER_HEIGHT_PX}px;
    right: 0px;
    z-index: 1;
  }
`;
