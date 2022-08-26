import styled from "@emotion/styled";
import { GRAY, PRIMARY_BLUE } from "src/components/common/theme";

export const CollectionMetadata = styled.div`
  display: grid;
  gap: 8px 24px;
  grid-area: metadata;
  grid-template-columns: 96px auto;
  justify-self: flex-start;
  letter-spacing: -0.1px;
  line-height: 18px;
`;

export const MetadataLabel = styled.span`
  color: ${GRAY.A};
`;

export const MetadataValue = styled.a`
  color: ${PRIMARY_BLUE};
  word-break: break-word;

  &:focus {
    outline: none;
  }

  &:hover {
    color: ${PRIMARY_BLUE};
  }
`;
