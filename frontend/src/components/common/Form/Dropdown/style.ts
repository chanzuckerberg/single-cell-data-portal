import styled from "@emotion/styled";
import { CommonThemeProps } from "@czi-sds/components";
import { gray400 } from "src/common/theme";

interface DropdownFormProps extends CommonThemeProps {
  isSelected?: boolean;
}

export const DropdownForm = styled("div")<DropdownFormProps>`
  grid-column: 1 / -1;

  .MuiButton-root {
    border-color: ${gray400};
    width: 100%;
  }
`;
