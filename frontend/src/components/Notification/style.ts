import styled from "@emotion/styled";
import {
  CommonThemeProps,
  Icon,
  fontBodyXs,
  fontHeaderXs,
  getColors,
  getCorners,
  getSpaces,
} from "@czi-sds/components";
import { spacesM, spacesXl } from "src/common/theme";

interface StyledNotificationDetailsProps extends CommonThemeProps {
  isCitation: boolean;
}

export const StyledNotificationWrapper = styled.div`
  position: absolute;
  z-index: 999;
  overflow: hidden;
`;

export const StyledIcon = styled(Icon)`
  height: 20px;
  width: 20px;
`;

export const StyledNotificationLabel = styled.div`
  ${fontHeaderXs}
  margin: 0 !important;
  padding-bottom: ${spacesM}px;
  color: black;
`;

export const StyledNotificationDetails = styled.div<StyledNotificationDetailsProps>`
  ${fontBodyXs}
  margin: 0 !important;
  color: black;
  width: 340px;
  height: auto;
  ${(props) => {
    const { isCitation } = props;
    const colors = getColors(props);
    const spaces = getSpaces(props);
    const corners = getCorners(props);
    const border = isCitation && `border: 1px solid ${colors?.gray[200]}`;
    const backgroundColor =
      isCitation && `background-color: ${colors?.gray[100]}`;
    const padding = isCitation && `padding: ${spaces?.m}px`;
    const borderRadius = isCitation && `border-radius: ${corners?.m}px`;
    return `
      ${border};
      ${backgroundColor};
      ${padding};
      ${borderRadius};
    `;
  }}
`;
