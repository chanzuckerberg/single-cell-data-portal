import styled from "@emotion/styled";
import { fontBodyXxs, getColors } from "@czi-sds/components";

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
