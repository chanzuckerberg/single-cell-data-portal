import { GetServerSideProps } from "next";
import { getServerSideSitemap, ISitemapField } from "next-sitemap";

// github copilot convert timestamp function
function convertTimestamp(timestamp: string) {
  const date = new Date(parseFloat(timestamp) * 1000);
  const year = date.getFullYear();
  const month = date.getMonth() + 1;
  const day = date.getDate();
  const hour = date.getHours();
  const minute = date.getMinutes();
  const second = date.getSeconds();
  const offset = date.getTimezoneOffset();
  const offsetHours = Math.abs(Math.floor(offset / 60));
  const offsetMinutes = Math.abs(offset % 60);
  const timezoneOffset = `${offset < 0 ? "+" : "-"}${offsetHours
    .toString()
    .padStart(2, "0")}:${offsetMinutes.toString().padStart(2, "0")}`;

  return `${year}-${month.toString().padStart(2, "0")}-${day
    .toString()
    .padStart(2, "0")}T${hour.toString().padStart(2, "0")}:${minute
    .toString()
    .padStart(2, "0")}:${second.toString().padStart(2, "0")}${timezoneOffset}`;
}

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
