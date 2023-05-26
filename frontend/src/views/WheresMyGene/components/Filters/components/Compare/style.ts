import styled from "@emotion/styled";
import { CommonThemeProps, getSpaces, Chip } from "@czi-sds/components";
import { PINK } from "src/common/theme";

export const NewChip = styled(Chip)`
  background-color: ${PINK};
  color: white;
  height: 20px !important;

  .MuiChip-label {
    ${(props: CommonThemeProps) => {
      const spaces = getSpaces(props);

      return `
        padding: ${spaces?.xxs}px ${spaces?.xs}px;
        font-weight: 500;
    `;
    }}
  }
`;

export const LabelWrapper = styled("div")`
  display: flex;

  ${(props: CommonThemeProps) => {
    const spaces = getSpaces(props);

    return `
      gap: ${spaces?.xs}px;
    `;
  }}
`;
