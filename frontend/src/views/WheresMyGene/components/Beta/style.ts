import styled from "@emotion/styled";
import { Callout, getColors } from "czifui";

export const StyledCallout = styled(Callout)`
  width: 650px;
  margin-top: 0;

  ${(props) => {
    const colors = getColors(props);

    return `
      background: ${colors?.beta[100]};

      .MuiAlert-icon {
        path {
          fill: ${colors?.beta[600]};
        }
      }
    `;
  }}
`;

export const SubmitIssue = styled.a`
  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.beta[600]};
    `;
  }}
`;
