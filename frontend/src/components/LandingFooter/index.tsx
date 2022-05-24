import Image from "next/image";
import Link from "next/link";
import { ROUTES } from "src/common/constants/routes";
import CZILogo from "src/components/common/staticPages/czi-logo-white.png";
import Logo from "src/components/common/staticPages/NEWxLOGO.png";
import styles from "./index.module.scss";

const LandingFooter = (): JSX.Element => {
  return (
    <footer className={styles.footer}>
      <div className={styles.footerTopContainer}>
        <div className={styles.footerLogo}>
          <Link href={ROUTES.HOMEPAGE} passHref>
            <a>
              <Image
                data-test-id="cellxgene-logo"
                src={Logo}
                alt="cellxgene logo"
              />
            </a>
          </Link>
        </div>

        <div className={styles.footerTopLinks}>
          <a
            href="https://github.com/chanzuckerberg/aspen/"
            target="_blank"
            rel="noopener"
          >
            Github
          </a>
          <a
            href="https://chanzuckerberg.com/careers/career-opportunities/?team=data,design,engineering,it,infrastructure-security,product,technical-program-management/"
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
