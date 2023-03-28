import styled from "@emotion/styled";
import { ButtonIcon, fontBodyXs, fontBodyXxs, getColors, Icon } from "czifui";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";

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

export const StyledNotificationWrapper = styled.div`
  position: absolute;
  top: ${HEADER_HEIGHT_PX}px;
  right: 24px;
  z-index: 999;
`;

export const StyledIcon = styled(Icon)`
  height: 20px;
  width: 20px;
`;

export const StyledNotificationLabel = styled.div`
  ${fontBodyXs}

  margin: 0 !important;

  color: black;
`;

export const StyledNotificationDetails = styled.div`
  ${fontBodyXs}

  margin: 0 !important;

  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.gray[500]}
    `;
  }}
`;
