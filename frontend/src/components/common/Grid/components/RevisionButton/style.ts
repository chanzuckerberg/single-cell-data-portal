import styled from "@emotion/styled";
import { ButtonIcon } from "@czi-sds/components";
import { primary600 } from "src/common/theme";

export const RevisionButton = styled(ButtonIcon)`
  &:hover {
    color: ${primary600};
  }

  &:active {
    color: ${primary600};
  }

  .MuiSvgIcon-root {
    height: 16px;
    width: 16px;
  }
`;
