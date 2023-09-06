import { GetServerSideProps } from "next";
import { getServerSideSitemap, ISitemapField } from "next-sitemap";
import { convertTimestamp } from "src/common/sitemaps/utils";
import { CELLGUIDE_DATA_URL } from "src/configs/configs";
import xml2js from "xml2js";

interface CellGuideData {
  ListBucketResult: {
    Contents: {
      Key: string;
      LastModified: string;
    };
  };
}

export const getServerSideProps: GetServerSideProps = async (context) => {
  const response = await fetch(CELLGUIDE_DATA_URL);

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
