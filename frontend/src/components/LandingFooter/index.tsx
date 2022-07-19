import Image from "next/image";
import { ROUTES } from "src/common/constants/routes";
import CZILogo from "src/components/common/staticPages/czi-logo-white.png";
import { HomepageLink } from "../common/HomepageLink";
import styles from "./index.module.scss";

const LandingFooter = (): JSX.Element => {
  return (
    <footer className={styles.footer}>
      <div className={styles.footerTopContainer}>
        <div className={styles.footerLogo}>
          <HomepageLink />
        </div>

        <div className={styles.footerTopLinks}>
          <a
            href="https://github.com/chanzuckerberg/single-cell-data-portal"
            target="_blank"
            rel="noopener"
          >
            Github
          </a>
          <a
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
          <a href="mailto:cellxgene@chanzuckerberg.com">Contact Us</a>
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
  );
};

export default LandingFooter;
