import Image from "next/image";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { ROUTES } from "src/common/constants/routes";
import wordmark from "src/common/images/cellxgene-discover-wordmark.svg";
import CZILogo from "src/components/common/staticPages/czi-logo-white.png";
import BottomBanner from "../BottomBanner";
import styles from "./index.module.scss";
import { GENE_EXPRESSION_BANNER_SURVEY_LINK } from "src/common/constants/airtableLinks";

const LandingFooter = (): JSX.Element => {
  return (
    <>
      <BottomBanner
        asFooter
        airtableLink={GENE_EXPRESSION_BANNER_SURVEY_LINK}
      />
      <footer className={styles.footer}>
        <div className={styles.footerTopContainer}>
          <div className={styles.footerLogo}>
            <a href={ROUTES.HOMEPAGE} target="_blank" rel="noopener">
              <Image src={wordmark} alt="CZ CELLxGENE Discover" />
            </a>
          </div>

          <div className={styles.footerTopLinks}>
            <a
              onClick={() => {
                track(EVENTS.GITHUB_CLICKED);
              }}
              href="https://github.com/chanzuckerberg/single-cell-data-portal"
              target="_blank"
              rel="noopener"
            >
              Github
            </a>
            <a
              onClick={() => {
                track(EVENTS.BROWSE_CAREERS_CLICKED, {
                  button: "careers",
                });
              }}
              href="https://chanzuckerberg.com/careers/career-opportunities/?team=data,design,engineering,product,technical-program-management&initiative=science&gh_src=20d9f28d1us"
              target="_blank"
              rel="noopener"
            >
              Careers
            </a>
          </div>
        </div>

        <div className={styles.footerBottomContainer}>
          <div className={styles.footerBottomLinks}>
            <a href={ROUTES.PRIVACY} target="_blank" rel="noopener">
              Privacy
            </a>
            <a href={ROUTES.TOS} target="_blank" rel="noopener">
              Terms
            </a>
            <a href={ROUTES.SITEMAP}>Sitemap</a>
            <a
              onClick={() => {
                track(EVENTS.CONTACT_US_CLICKED);
              }}
              href="mailto:cellxgene@chanzuckerberg.com"
            >
              Contact Us
            </a>
          </div>

          <div className={styles.footerBottomLogos}>
            <span className={styles.footerBottomLogoText}>
              In partnership with:
            </span>
            <div className={styles.footerBottomLogosInner}>
              <a
                href="https://chanzuckerberg.com/"
                target="_blank"
                rel="noopener"
                className={styles.footerBottomLogoLeft}
              >
                <Image src={CZILogo} alt="Chan Zuckerberg Initiative logo" />
              </a>
            </div>
          </div>
        </div>
      </footer>
    </>
  );
};

export default LandingFooter;
