import loadable from "@loadable/component";
import Head from "next/head";
import { FC } from "react";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import Collections from "src/components/Collections";

const AsyncUploadCSV = loadable(
  () =>
    /*webpackChunkName: 'src/components/UploadCSV' */ import(
      "src/components/UploadCSV"
    )
);
const Homepage: FC = () => {
  // (thuang): TEMP. Remove when we do https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/corpora-data-portal/917
  const isGeneSetsOn = get(FEATURES.GENE_SETS) === BOOLEAN.TRUE;

  return (
    <>
      <Head>
        <title>cellxgene | Homepage</title>
      </Head>
      {
        // (thuang): TEMP. Remove when we do https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/corpora-data-portal/917
        isGeneSetsOn && <AsyncUploadCSV />
      }
      <Collections />
    </>
  );
};

export default Homepage;
