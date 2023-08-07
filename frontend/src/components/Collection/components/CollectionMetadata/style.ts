import styled from "@emotion/styled";
import { fontBodyS } from "@czi-sds/components";
import { gray500, primary400, spacesS, spacesXl } from "src/common/theme";

export const CollectionMetadata = styled.div`
  display: grid;
  ${fontBodyS}
  gap: ${spacesS}px ${spacesXl}px;
  grid-area: metadata;
  grid-template-columns: 96px auto;
  justify-self: flex-start;
  letter-spacing: -0.006em;
`;

export const MetadataLabel = styled.span`
  color: ${gray500};
`;

export const MetadataValue = styled.a`
  color: ${primary400};
  word-break: break-word;

  &:focus {
    outline: none;
  }

  &:hover {
    color: ${primary400};
  }
`;
