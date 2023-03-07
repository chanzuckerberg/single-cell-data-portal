import styled from "@emotion/styled";
import { CommonThemeProps, getSpaces, Tag } from "czifui";
import { PINK } from "src/common/theme";

export const NewTag = styled(Tag)`
  background-color: ${PINK};
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
