import { FC } from "react";
import {
  INCOMPLETE_GENE_MIGRATIONS_COLLECTIONS,
  INCOMPLETE_GENE_MIGRATIONS_HEADER,
  INCOMPLETE_GENE_MIGRATIONS_MSG,
} from "./constants";
import { BannerHeader, BannerWrapper, Wrapper } from "./style";

interface Props {
  collection_id: string;
}

const Banner: FC<Props> = ({ collection_id }) => {
  if (!INCOMPLETE_GENE_MIGRATIONS_COLLECTIONS.includes(collection_id))
    return null;

  const bannerHeader = INCOMPLETE_GENE_MIGRATIONS_HEADER;
  const bannerMsg = INCOMPLETE_GENE_MIGRATIONS_MSG;

  return (
    <Wrapper>
      <BannerHeader>{bannerHeader}</BannerHeader>
      <BannerWrapper>{bannerMsg}</BannerWrapper>
    </Wrapper>
  );
};

export default Banner;
