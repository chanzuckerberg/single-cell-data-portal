import { GetServerSideProps } from "next";
import { getServerSideSitemap, ISitemapField } from "next-sitemap";
import { convertTimestamp } from "src/common/sitemaps/utils";

export const getServerSideProps: GetServerSideProps = async (context) => {
  const response = await fetch(
    "https://api.cellxgene.cziscience.com/dp/v1/collections/index"
  );

  const collections: any[] = await response.json();

  const fields: ISitemapField[] = collections.map((collection) => ({
    loc: `https://www.cellxgene.cziscience.com/collections/${collection.id}`,
    lastmod: collection.revised_at
      ? convertTimestamp(collection.revised_at)
      : convertTimestamp(collection.published_at),
  }));

  return getServerSideSitemap(context, fields);
};

export default function Site() {
  // intentionally empty, since we're generating xml files
}
