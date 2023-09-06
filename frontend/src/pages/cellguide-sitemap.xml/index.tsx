import { GetServerSideProps } from "next";
import { getServerSideSitemap, ISitemapField } from "next-sitemap";
import { convertTimestamp } from "src/common/sitemaps/utils";
import { CELLGUIDE_DATA_URL } from "src/configs/configs";

export const getServerSideProps: GetServerSideProps = async (context) => {
  const latestSnapshotRes = await fetch(
    `${CELLGUIDE_DATA_URL}/latest_snapshot_identifier`
  );

  const latestSnapshotId = await latestSnapshotRes.json();

  const response = await fetch(
    `${CELLGUIDE_DATA_URL}/${latestSnapshotId}/celltype_metadata.json`
  );

  if (!response.ok) {
    throw new Error(`Failed to fetch data. Status: ${response.status}`);
  }

  const cellGuideData = await response.json();

  const cellGuideIds = Object.keys(cellGuideData);

  const fields: ISitemapField[] = cellGuideIds.map((id: string) => {
    const date = new Date().getTime();

    return {
      loc: `https://cellxgene.cziscience.com/cellguide/${id}`,
      lastmod: convertTimestamp(date.toString()),
    };
  });

  return getServerSideSitemap(context, fields);
};

export default function Site() {
  // intentionally empty, since we're generating xml files
}
