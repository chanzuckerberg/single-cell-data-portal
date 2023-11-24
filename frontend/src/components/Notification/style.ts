import styled from "@emotion/styled";
import {
  CommonThemeProps,
  Icon,
  fontBodyXs,
  fontHeaderXs,
} from "@czi-sds/components";
import { NotificationWrapper } from "src/components/common/Filter/common/style";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";

interface StyledNotificationDetailsProps extends CommonThemeProps {
  isCitation: boolean;
}

export const StyledNotificationWrapper = styled(NotificationWrapper)`
  top: ${HEADER_HEIGHT_PX}px;
  right: 24px;
`;

export const StyledIcon = styled(Icon)`
  height: 20px;
  width: 20px;
`;

export const StyledNotificationLabel = styled.div`
  ${fontHeaderXs}
  margin: 0 !important;
  color: black;
`;

export const StyledNotificationDetails = styled.div<StyledNotificationDetailsProps>`
  ${fontBodyXs}
  margin: 0 !important;
  color: black;
  width: 340px;
  height: auto;
  padding: 8px 12px 8px 12px;
  border-radius: 8px;
  ${(props) => {
    const { isCitation } = props;
    const border = isCitation ? ` 1px solid #EAEAEA` : null;
    const backgroundColor = isCitation ? ` #F8F8F8` : null;
    console.log("isCitation", isCitation);
    console.log("border", border);
    console.log("backgroundColor", backgroundColor);

    return `
      border: ${border};
      background-color: ${backgroundColor};
    `;
  }}
`;
