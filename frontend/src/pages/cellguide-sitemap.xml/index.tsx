import { GetServerSideProps } from "next";
import { getServerSideSitemap, ISitemapField } from "next-sitemap";
import { convertTimestamp } from "src/common/sitemaps/utils";
import { CELLGUIDE_DATA_URL } from "src/configs/configs";
import { ROUTES } from "src/common/constants/routes";

export const getServerSideProps: GetServerSideProps = async (context) => {
  const latestSnapshotRes = await fetch(
    `${CELLGUIDE_DATA_URL}/latest_snapshot_identifier`
  );

  const latestSnapshotId = await latestSnapshotRes.json();

  const response = await fetch(
    `${CELLGUIDE_DATA_URL}/${latestSnapshotId}/celltype_metadata.json`
  );

  const tissueResponse = await fetch(
    `${CELLGUIDE_DATA_URL}/${latestSnapshotId}/tissue_metadata.json`
  );

  if (!response.ok) {
    throw new Error(`Failed to fetch data. Status: ${response.status}`);
  }

  if (!tissueResponse.ok) {
    throw new Error(`Failed to fetch data. Status: ${tissueResponse.status}`);
  }

  const cellGuideData = await response.json();

  const tissueData = await tissueResponse.json();

  const cellGuideIds = Object.keys(cellGuideData);

  const tissueDataIds = Object.keys(tissueData);

  const allIds = [...cellGuideIds, ...tissueDataIds];

  const fields: ISitemapField[] = allIds.map((id: string) => {
    const date = Math.floor(new Date().getTime() / 1000);

    const endpoint = id.includes("UBERON")
      ? `https://cellxgene.cziscience.com${ROUTES.CELL_GUIDE}/tissues/${id}`
      : `https://cellxgene.cziscience.com${ROUTES.CELL_GUIDE}/${id}`;

    return {
      loc: endpoint,
      lastmod: convertTimestamp(date.toString()),
    };
  });

  return getServerSideSitemap(context, fields);
};

export default function Site() {
  // intentionally empty, since we're generating xml files
}
