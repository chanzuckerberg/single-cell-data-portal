import {
  BannerHeader,
  BannerText,
  CollectionBanner,
  Wrapper,
} from "src/views/Collection/components/Banner/style";
import {
  INCOMPLETE_GENE_MIGRATIONS_COLLECTIONS,
  INCOMPLETE_GENE_MIGRATIONS_HEADER,
  INCOMPLETE_GENE_MIGRATIONS_MSG,
} from "./constants";

interface Props {
  collectionId: string;
  isFilterEnabled?: boolean;
}

const Banner = ({
  collectionId,
  isFilterEnabled = false,
}: Props): JSX.Element | null => {
  if (!INCOMPLETE_GENE_MIGRATIONS_COLLECTIONS.includes(collectionId))
    return null;

  const BannerWrapper = isFilterEnabled ? CollectionBanner : Wrapper;
  const bannerHeader = INCOMPLETE_GENE_MIGRATIONS_HEADER;
  const bannerMsg = INCOMPLETE_GENE_MIGRATIONS_MSG;

  return (
    <BannerWrapper>
      <BannerHeader>{bannerHeader}</BannerHeader>
      <BannerText>{bannerMsg}</BannerText>
    </BannerWrapper>
  );
};

export default Banner;
