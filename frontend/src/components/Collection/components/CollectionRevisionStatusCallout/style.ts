import styled from "@emotion/styled";
import { Banner, Link } from "@czi-sds/components";
import { CollectionHero } from "src/views/Collection/style";

export const CollectionRevisionCallout = styled(Banner)`
  height: 39px;
  left: 0;
  letter-spacing: -0.006em;
  position: absolute;
  top: 0;

  > div > div {
    display: contents; /* targeting @czi-sds/components IconWrapper style to override iconSize "l" and spaces "m" setting in app theme */
  }

  svg {
    height: 22px;
    margin-right: 10px;
    width: 22px;
  }

  + ${CollectionHero} {
    /* CollectionView top padding is 40px and CollectionRevisionCallout height is 39px.
      Therefore to meet the desired spacing of 24px between CollectionRevisionCallout and CollectionHero,
      margin-top should only be 23px. */
    margin-top: 23px;
  }
`;

export const TextLink = styled(Link)`
  color: inherit;
  font-weight: 500;
  text-decoration: underline;

  &:hover {
    color: inherit;
  }
`;
