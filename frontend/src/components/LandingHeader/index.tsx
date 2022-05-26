import Image from "next/image";
import Link from "next/link";
import { useRouter } from "next/router";
import { ROUTES } from "src/common/constants/routes";
import Logo from "src/components/common/staticPages/NEWxLOGO.png";
import styles from "./index.module.scss";

const LandingHeader = (): JSX.Element => {
  const router = useRouter();

  return (
    <header className={styles.header}>
      <Link href={ROUTES.HOMEPAGE} passHref>
        <a className={styles.headerLogoContainer}>
          <Image
            data-test-id="cellxgene-logo"
            src={Logo}
            alt="cellxgene logo"
          />
        </a>
      </Link>
      <nav className={styles.headerNavContainer}>
        <div className={styles.headerNavLeft}>
          {/* href update tbd */}
          <Link href={ROUTES.HOMEPAGE} passHref>
            <a>Browse</a>
          </Link>
          {/* href update tbd */}
          <Link href={ROUTES.HOMEPAGE} passHref>
            <a>Atlas</a>
          </Link>
          {/* href update tbd */}
          <Link href={ROUTES.HOMEPAGE} passHref>
            <a>Explore</a>
          </Link>
        </div>
        <div className={styles.headerNavLeft}>
          {/* CHANGE TO router.route === "/" FOR PROD */}
          {router.route === "/landing-page" && (
            <Link
              href="https://chanzuckerberg.com/careers/career-opportunities/?team=data,design,engineering,product,technical-program-management&initiative=science&gh_src=20d9f28d1us"
              passHref
            >
              <a className={styles.btnLink}>We're Hiring</a>
            </Link>
          )}
          {/* href update tbd */}
          <Link href={ROUTES.HOMEPAGE} passHref>
            <a>Help & Documentation</a>
          </Link>
          {/* href update tbd */}
          <Link href={ROUTES.HOMEPAGE} passHref>
            <a>Login</a>
          </Link>
        </div>
      </nav>
    </header>
  );
};

export default LandingHeader;
