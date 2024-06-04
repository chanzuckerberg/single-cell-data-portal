import styled from "@emotion/styled";
import { View } from "src/views/globalStyle";
import { spacesXl } from "src/common/theme";

export const DatasetsView = styled(View)`
  display: grid;
  gap: ${spacesXl}px;
  place-content: flex-start stretch;
`;
