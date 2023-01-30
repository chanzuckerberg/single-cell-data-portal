import { GetServerSideProps } from "next";
import { getServerSideSitemap, ISitemapField } from "next-sitemap";

// function to convert epoch/unix timestamp to YYYY-MM-DD HH:MM +00:00 format
function unixToYYYYMMDD(date: string) {
  let myDate = new Date(parseFloat(date) * 1000);

  let month = String(myDate.getMonth() + 1);
  if (parseFloat(month) < 10) {
    month = "0" + String(month);
  }

  let day = String(myDate.getDate());
  if (parseFloat(day) < 10) {
    day = "0" + String(day);
  }

  let hours = String(myDate.getHours());
  if (parseFloat(hours) < 10) {
    hours = "0" + String(hours);
  }

  let minutes = String(myDate.getMinutes());
  if (parseFloat(minutes) < 10) {
    minutes = "0" + String(minutes);
  }

  let dateStr =
    myDate.getFullYear() +
    "-" +
    month +
    "-" +
    day +
    " " +
    hours +
    ":" +
    minutes +
    " +00:00";

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
