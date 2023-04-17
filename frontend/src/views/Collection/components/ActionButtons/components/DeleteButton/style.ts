import { CommonThemeProps, getColors } from "czifui";
import styled from "@emotion/styled";
import { ActionButton } from "src/views/Collection/components/ActionButtons/style";

const error400 = (props: CommonThemeProps) => getColors(props)?.error[400];
const error500 = (props: CommonThemeProps) => getColors(props)?.error[500];

export const DeleteButton = styled(ActionButton)`
  background-color: ${error400};

  &:hover {
    background-color: ${error500};
  }
`;
