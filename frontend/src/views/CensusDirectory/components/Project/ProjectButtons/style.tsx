import { Button } from "@czi-sds/components";
import styled from "@emotion/styled";
import { fontWeightMedium, spacesDefault, spacesXl } from "src/common/theme";

export const ButtonsColumn = styled.div`
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
`;

export const ButtonsRow = styled.div`
  display: flex;
  flex-direction: row;
  gap: ${spacesDefault}px;
  margin-bottom: ${spacesXl}px;

  // item container
  & > div {
    margin-right: ${spacesXl}px;
  }
`;

export const StyledButton = styled(Button)`
  font-weight: ${fontWeightMedium};
  min-width: 80px;
`;
