import React, { FC } from "react";
import rawCellxgeneLogo from "src/components/common/staticPages/cellxgene.png";
import rawCZILogo from "src/components/common/staticPages/CZI_Logotype_RGB.png";
import {
  CellxgeneLogo,
  CommonStyle,
  CZILogo,
  Layout,
  PrivacyStyle,
} from "src/components/common/staticPages/style";
import SEO from "src/components/seo";

const PreviewPolicies: FC = () => {
  return (
    <Layout>
      <CommonStyle>
        <PrivacyStyle>
          <SEO title="Preview Policies" />
          <header>
            <CZILogo
              data-test-id="czi-logo"
              src={String(rawCZILogo)}
              alt="CZI logo"
            />
            <CellxgeneLogo
              data-test-id="cellxgene-logo"
              src={String(rawCellxgeneLogo)}
              alt="cellxgene logo"
            />
          </header>

          <h1>
            Summary of April 2021 Updates to Terms of Service and Privacy Policy
          </h1>

          <p>
            <span>
              We periodically update our terms and policies to ensure that they
              are transparent and easy to understand. Here is a summary of the
              key updates from April 2021:
              <ul>
                <li>
                  Cellxgene is now provided by the Chan Zuckerberg Initiative
                  Foundation. The Chan Zuckerberg Initiative has moved its
                  Science team to its 501(c)(3) private foundation from CZI,
                  LLC. We’re committed to preserving the same experience for
                  you. For questions about these changes, please contact us at{" "}
                  <a href="mailto:privacy@chanzuckerberg.com">
                    privacy@chanzuckerberg.com
                  </a>
                  .
                </li>
              </ul>
            </span>
          </p>
        </PrivacyStyle>
      </CommonStyle>
    </Layout>
  );
};

export default PreviewPolicies;
