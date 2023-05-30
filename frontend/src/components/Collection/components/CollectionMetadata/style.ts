import styled from "@emotion/styled";
import {
  CommonThemeProps,
  fontBodyS,
  getColors,
  getSpaces,
} from "@czi-sds/components";

const gray500 = (props: CommonThemeProps) => getColors(props)?.gray[500];
const primary400 = (props: CommonThemeProps) => getColors(props)?.primary[400];
const spacesS = (props: CommonThemeProps) => getSpaces(props)?.s;
const spacesXl = (props: CommonThemeProps) => getSpaces(props)?.xl;

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
