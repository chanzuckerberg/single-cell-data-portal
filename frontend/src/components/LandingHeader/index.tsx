import Image from "next/image";
import Link from "next/link";
import { useRouter } from "next/router";
import Logo from "src/components/common/staticPages/NEWxLOGO.png";
import styles from "./index.module.scss";

const LandingHeader = (): JSX.Element => {
  const router = useRouter();

  console.log(router);

  return (
    <header className={styles.header}>
      <Link href="/" passHref>
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
          <Link href="/" passHref>
            <a>Browse</a>
          </Link>
          <Link href="/" passHref>
            <a>Atlas</a>
          </Link>
          <Link href="/" passHref>
            <a>Explore</a>
          </Link>
        </div>
        <div className={styles.headerNavLeft}>
          {/* CHANGE TO router.route === "/" FOR PROD */}
          {router.route === "/landing-page" && (
            <Link
              href="https://chanzuckerberg.com/careers/career-opportunities/?team=data,design,engineering,it,infrastructure-security,product,technical-program-management/"
              passHref
            >
              <a className={styles.btnLink}>We're Hiring</a>
            </Link>
          )}
          <Link href="/" passHref>
            <a>Help & Documentation</a>
          </Link>
          <Link href="/" passHref>
            <a>Login</a>
          </Link>
        </div>
      </nav>
    </header>
  );
};

export default LandingHeader;
