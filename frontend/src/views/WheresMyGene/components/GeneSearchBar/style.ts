import styled from "@emotion/styled";
import { Button, fontBodyXxs, getColors } from "@czi-sds/components";

export const Container = styled.div`
  width: 80vw;
`;

export const ActionWrapper = styled.div`
  display: flex;
  gap: 16px;
`;

export const Label = styled.label`
  ${fontBodyXxs}

  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.gray[500]}
    `;
  }}
`;

export const LoadingIndicatorWrapper = styled.div`
  display: flex;
  align-items: center;
`;

export const StyledButtonWrapper = styled.div`
  align-self: center;
`;

export const StyledClearButton = styled(Button)`
  white-space: nowrap;
  font-weight: 500;
`;
