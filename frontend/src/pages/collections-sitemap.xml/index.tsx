import { GetServerSideProps } from "next";
import { getServerSideSitemap, ISitemapField } from "next-sitemap";

// function to convert epoch/unix timestamp to YYYY-MM-DD format
function unixToYYYYMMDD(date: string) {
  let myDate = new Date(parseFloat(date) * 1000);

  let dateStr =
    myDate.getFullYear() +
    "-" +
    (myDate.getMonth() + 1) +
    "-" +
    myDate.getDate();

  return dateStr;
}

export const getServerSideProps: GetServerSideProps = async (context) => {
  const response = await fetch(
    "https://api.cellxgene.cziscience.com/dp/v1/collections/index"
  );

  const collections: any[] = await response.json();

  const fields: ISitemapField[] = collections.map((collection) => ({
    loc: `https://www.cellxgene.cziscience.com/collections/${collection.id}`,
    lastmod: collection.revised_at
      ? unixToYYYYMMDD(collection.revised_at)
      : unixToYYYYMMDD(collection.published_at),
  }));

  return getServerSideSitemap(context, fields);
};

export default function Site() {}
