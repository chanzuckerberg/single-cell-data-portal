import styled from "@emotion/styled";
import { fontBodyXs, Tag as SDSTag } from "@czi-sds/components";
import {
  primary100,
  primary500,
  spacesS,
  spacesXxs,
  spacesXxxs,
  success100,
  success500,
} from "src/common/theme";

export const StatusTags = styled.div`
  display: grid;
  gap: ${spacesXxs}px;
  justify-items: flex-start;
`;

export const StatusTag = styled.div`
  align-items: center;
  display: flex;
  gap: ${spacesXxs}px;
`;

export const Tag = styled(SDSTag)`
  &.MuiChip-root {
    cursor: default;
    margin: 0;
    padding: ${spacesXxxs}px ${spacesS}px;
  }

  .MuiChip-label {
    ${fontBodyXs}
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
