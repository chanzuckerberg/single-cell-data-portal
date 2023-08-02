import styled from "@emotion/styled";
import { ButtonIcon, fontBodyXs, fontBodyXxs, Icon } from "@czi-sds/components";
import { NotificationWrapper } from "src/components/common/Filter/common/style";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";
import { gray500 } from "src/common/theme";

export const StyledButtonDiv = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
  padding: 0 10px;
`;

export const StyledLabel = styled.div`
  ${fontBodyXxs}

  white-space: nowrap;

  color: ${gray500};
`;

export const StyledButtonIcon = styled(ButtonIcon)`
  width: 30px;
  height: 30px;
`;

export const StyledNotificationWrapper = styled(NotificationWrapper)`
  top: ${HEADER_HEIGHT_PX}px;
  right: 24px;
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
  color: ${gray500};
`;
