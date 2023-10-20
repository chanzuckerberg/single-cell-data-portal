import styled from "@emotion/styled";
import { Chip } from "@czi-sds/components";
import { PINK, spacesXs, spacesXxs } from "src/common/theme";

export const NewChip = styled(Chip)`
  background-color: ${PINK};
  color: white;
  height: 20px !important;

  .MuiChip-label {
    padding: ${spacesXxs}px ${spacesXs}px;
    font-weight: 500;
  }
`;

export const LabelWrapper = styled("div")`
  display: flex;
  gap: ${spacesXs}px;
`;
