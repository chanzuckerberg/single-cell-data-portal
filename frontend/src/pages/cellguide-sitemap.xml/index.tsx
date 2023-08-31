import { GetServerSideProps } from "next";
import { getServerSideSitemap, ISitemapField } from "next-sitemap";
import xml2js from "xml2js";

function convertTimestamp(originalTimestamp: string) {
  const date = new Date(originalTimestamp);
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

interface CellGuideData {
  ListBucketResult: {
    Contents: {
      Key: string;
      LastModified: string;
    };
  };
}

export const getServerSideProps: GetServerSideProps = async (context) => {
  const response = await fetch(
    "https://cellguide.cellxgene.dev.single-cell.czi.technology"
  );

  const xmlData = await response.text();

  const cellTypeData = new Promise<CellGuideData>((resolve, reject) => {
    xml2js.parseString(xmlData, (err, jsonData) => {
      if (err) {
        console.error("Error parsing XML:", err);
        reject(err);
      } else {
        resolve(jsonData);
      }
    });
  });

  const cellTypeEntriesArray = Object.values(
    (await cellTypeData).ListBucketResult.Contents
  );

  const fields: ISitemapField[] = cellTypeEntriesArray.map((entry: any) => {
    const endpointSplitArr = entry.Key.toString().split("/");

    const endpointDotJson = endpointSplitArr[endpointSplitArr.length - 1];

    const endpoint = endpointDotJson.split(".json")[0];

    return {
      loc: `https://www.cellxgene.cziscience.com/cellguide/${endpoint}`,
      lastmod: convertTimestamp(entry.LastModified),
    };
  });

  return getServerSideSitemap(context, fields);
};

export default function Site() {
  // intentionally empty, since we're generating xml files
}
