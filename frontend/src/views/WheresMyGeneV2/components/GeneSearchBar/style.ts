import styled from "@emotion/styled";
import { fontBodyXxs } from "@czi-sds/components";
import { gray500 } from "src/common/theme";

export const Container = styled.div`
  width: 80vw;
`;

export const ActionWrapper = styled.div`
  display: flex;
  gap: 16px;
`;

export const Label = styled.label`
  ${fontBodyXxs}

  color: ${gray500};
`;

export const LoadingIndicatorWrapper = styled.div`
  display: flex;
  align-items: center;
`;
