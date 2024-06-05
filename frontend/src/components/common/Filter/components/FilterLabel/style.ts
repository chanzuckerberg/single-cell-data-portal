import styled from "@emotion/styled";
import { InputDropdown as SDSInputDropdown } from "@czi-sds/components";
import { spacesXxs } from "src/common/theme";

export const InputDropdown = styled(SDSInputDropdown)`
  &.MuiButton-root {
    padding-bottom: ${spacesXxs}px;
    padding-top: ${spacesXxs}px;

    .styled-label {
      font-weight: 500;
    }
  }
`;
