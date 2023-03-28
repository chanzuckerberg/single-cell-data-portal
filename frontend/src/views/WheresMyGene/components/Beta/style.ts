import styled from "@emotion/styled";
import { CommonThemeProps, getColors, Notification } from "czifui";

export const StyledNotification = styled(Notification)`
  ${(props: CommonThemeProps) => {
    const colors = getColors(props);

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

export const SubmitIssue = styled.a`
  ${(props: CommonThemeProps) => {
    const colors = getColors(props);

    return `
      color: ${colors?.beta[600]};
    `;
  }}
`;
