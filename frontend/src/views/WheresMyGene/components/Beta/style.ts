import styled from "@emotion/styled";
import { CommonThemeProps, getColors, Notification } from "czifui";

export const StyledNotification = styled(Notification)`
  ${(props: CommonThemeProps) => {
    const colors = getColors(props);

    return `
    .elevated {
      border-color: ${colors?.beta[400]} !important;
    }

    .MuiAlert-root {
      background-color: ${colors?.beta[100]};

      .MuiAlert-icon {
        path {
          fill: ${colors?.beta[600]};
        }
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
