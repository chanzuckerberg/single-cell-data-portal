import styled from "@emotion/styled";
import { CommonThemeProps } from "@czi-sds/components";
import { gray400, gray500, textPrimary } from "src/common/theme";

interface DropdownFormProps extends CommonThemeProps {
  isSelected?: boolean;
}

export const DropdownForm = styled("div")<DropdownFormProps>`
  .MuiButton-root {
    border-color: ${gray400};

    :hover {
      border-color: ${gray400};
    }

    height: 32px;
    width: 100%;
  }

  .MuiButton-text {
    > span {
      letter-spacing: -0.003em;

      ${(props) =>
        props.isSelected &&
        `
          color: ${textPrimary(props)};
        `}
    }

    svg {
      height: 16px;
      width: 16px;

      path {
        ${(props) => (!props.isSelected ? `fill: ${gray500};` : "")}
      }
    }
  }
`;
