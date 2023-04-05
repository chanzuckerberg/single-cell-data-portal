import styled from "@emotion/styled";
import { CommonThemeProps, fontBodyXs, getColors, Tag as SDSTag } from "czifui";

const primary100 = (props: CommonThemeProps) => getColors(props)?.primary[100];
const primary500 = (props: CommonThemeProps) => getColors(props)?.primary[500];
const success100 = (props: CommonThemeProps) => getColors(props)?.success[100];
const success500 = (props: CommonThemeProps) => getColors(props)?.success[500];

export const StatusTags = styled.div`
  display: grid;
  gap: 4px;
  justify-items: flex-start;
`;

export const StatusTag = styled.div`
  align-items: center;
  display: flex;
  gap: 4px;
`;

export const Tag = styled(SDSTag)`
  &.MuiChip-root {
    cursor: default;
    margin: 0;
    padding: 2px 8px;
  }

  .MuiChip-label {
    ${fontBodyXs};
    letter-spacing: -0.003em;
    text-transform: capitalize;
  }
`;

export const PrivateTag = styled(Tag)`
  &.MuiChip-root {
    background-color: ${primary100};

    .MuiChip-label {
      color: ${primary500};
    }
  }
`;

export const PublishedTag = styled(Tag)`
  &.MuiChip-root {
    background-color: ${success100};

    .MuiChip-label {
      color: ${success500};
    }
  }
`;

export const RevisionTag = styled(Tag)`
  &.MuiChip-root {
    background-color: ${primary100};

    .MuiChip-label {
      color: ${primary500};
    }
  }
`;
