import styled from "@emotion/styled";
import { Tag, fontBodyXs } from "@czi-sds/components";

export const StyledTag = styled(Tag)`
  ${fontBodyXs}
  & .MuiChip-label {
    color: #000000;
  }
`;

export const SourceDataTableWrapper = styled.div`
  max-width: calc(100vw - 32px);
`;
