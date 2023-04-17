import styled from "@emotion/styled";
import { CommonThemeProps, fontBodyS, getColors } from "czifui";

const gray500 = (props: CommonThemeProps) => getColors(props)?.gray[500];
const primary400 = (props: CommonThemeProps) => getColors(props)?.primary[400];

export const CollectionMetadata = styled.div`
  display: grid;
  ${fontBodyS}
  gap: 8px 24px;
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
