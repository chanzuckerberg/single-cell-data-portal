import styled from "@emotion/styled";
import { CommonThemeProps, getColors, Notification } from "czifui";

export const StyledNotification = styled(Notification)`
  ${(props: CommonThemeProps) => {
    const colors = getColors(props);

    // beta intent does not exist for SDS Notification, but the colors do
    return `
      #beta-notification {
        border-color: ${colors?.beta[400]} !important;
        background-color: ${colors?.beta[100]};

        path {
          fill: ${colors?.beta[600]};
        }
      }
    `;
  }}
`;

export const StyledNotificationWrapper = styled.div`
  bottom: 10px;
  position: absolute;
  right: 30px;
  width: 360px;
  z-index: 99;
`;

export const SubmitIssue = styled.a`
  ${(props: CommonThemeProps) => {
    const colors = getColors(props);

    return `
      color: ${colors?.beta[600]};
    `;
  }}
`;
