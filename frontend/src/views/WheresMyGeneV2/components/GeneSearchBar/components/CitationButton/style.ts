import styled from "@emotion/styled";
import { ButtonIcon, fontBodyXxs } from "@czi-sds/components";
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
  top: 4px;
`;
