import styled from "@emotion/styled";
import { fontBodyS, fontHeaderXl } from "@czi-sds/components";
import { spacesXxl } from "src/common/theme";

export const CollectionView = styled.div`
  padding: ${spacesXxl}px;
  position: relative; /* positions banner */
`;

export const CollectionHero = styled.div`
  align-items: flex-start; /* top aligns collection name with action buttons */
  column-gap: ${spacesXxl}px;
  display: flex;

  h3 {
    ${fontHeaderXl}
    flex: 1;
    letter-spacing: -0.019em;
    margin: 0;
  }
`;

export const CollectionDetail = styled.div`
  align-items: flex-start;
  display: grid;
  grid-template-areas: "consortia . ." "description . metadata"; /* grid areas for collection consortia and collection description and metadata with a null cell token (unnamed area) for the allocation of a gutter between the two columns */
  grid-template-columns: 8fr 1fr 8fr; /* grid columns for collection description and metadata (with 1fr allocated to column separation) */
  margin: 16px 0 44px;
`;

export const CollectionConsortia = styled("div")`
  ${fontBodyS}
  font-weight: 500;
  grid-area: consortia;
  letter-spacing: -0.006em;
  margin: -8px 0 8px; /* CollectionDetail margin top reduced to 8px by CollectionConsortia negative margin top */
`;
