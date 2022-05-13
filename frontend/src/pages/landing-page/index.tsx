import Head from "next/head";
import Image from "next/image";
import Link from "next/link";
import { useState } from "react";
import { useInView } from "react-intersection-observer";
import AnalyzeDatasetsImg from "src/components/common/staticPages/analyze-datasets.png";
import BiohubLogo from "src/components/common/staticPages/biohub-logo-white.png";
import HeroBg from "src/components/common/staticPages/cellxgene_hero_bg.png";
import LaptopImg from "src/components/common/staticPages/cellxgene_laptop.png";
import CZILogo from "src/components/common/staticPages/czi-logo-white.png";
import NewsImage1 from "src/components/common/staticPages/czi-news-tweet-1.png";
import NewsImage2 from "src/components/common/staticPages/czi-news-tweet-2.png";
import DownloadDataImg from "src/components/common/staticPages/download-data.png";
import ExpediteCollaborationImg from "src/components/common/staticPages/expedite-collaboration.png";
import GeneExpressionImg from "src/components/common/staticPages/gene-expression.png";
import Logo from "src/components/common/staticPages/NEWxLOGO.png";
import SingleCellDataImg from "src/components/common/staticPages/single-cell-data.png";
import AnalyzeDatasetsIcon from "./icons/analyze-datasets";
import DownloadDataIcon from "./icons/download-data";
import ExpediteCollaborationIcon from "./icons/expedite-collaboration";
import LinkArrow from "./icons/external-link-arrow";
import GeneExpressionIcon from "./icons/gene-expression";
import SingleCellDataIcon from "./icons/single-cell-data";
import styles from "./index.module.scss";

const LandingPage = (): JSX.Element => {
  const { ref: section1, inView: inView1 } = useInView({
    rootMargin: "-70% 0px -30% 0px",
  });
  const { ref: section2, inView: inView2 } = useInView({
    rootMargin: "-70% 0px -30% 0px",
  });
  const { ref: section3, inView: inView3 } = useInView({
    rootMargin: "-70% 0px -30% 0px",
  });
  const { ref: section4, inView: inView4 } = useInView({
    rootMargin: "-70% 0px -30% 0px",
  });
  const { ref: section5, inView: inView5 } = useInView({
    rootMargin: "-70% 0px -30% 0px",
  });

  const [activeNav, setActiveNav] = useState(1);

  return (
    <>
      <Head>
        <title>cellxgene | Home</title>
      </Head>
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
            <Link href="/" passHref>
              <a className={styles.btnLink}>We're Hiring</a>
            </Link>
            <Link href="/" passHref>
              <a>Help & Documentation</a>
            </Link>
            <Link href="/" passHref>
              <a>Login</a>
            </Link>
          </div>
        </nav>
      </header>
      <div className={styles.heroContainer}>
        <div
          className={styles.heroImgContainer}
          style={{
            backgroundImage: `url(${HeroBg.src})`,
            backgroundSize: "cover",
            backgroundPosition: "center center",
            backgroundRepeat: "no-repeat",
          }}
        >
          <div className={styles.laptopImg}>
            <Image src={LaptopImg} alt="" />
          </div>
        </div>
        <div className={styles.heroTextContainer}>
          <h1>Discover the mechanisms of human health</h1>
          <p>
            Download and visually explore reference-quality data to understand
            the functionality of human tissues at the cellular level.
          </p>
          <div className={styles.heroStatsContainer}>
            <div>
              <span>Cells</span>
              <p>100M+</p>
            </div>
            <div>
              <span>datasets</span>
              <p>436</p>
            </div>
            <div>
              <span>donors</span>
              <p>2.7k+</p>
            </div>
          </div>
        </div>
      </div>
      <main className={styles.main}>
        <div className={styles.contentNav}>
          <Link passHref href="#single-cell">
            <a
              onClick={() => {
                setActiveNav(1);
              }}
              className={`${styles.contentLink} ${
                activeNav === 1 || inView1 ? styles.active : ""
              }`}
            >
              Find single-cell data
            </a>
          </Link>
          <Link passHref href="#gene-expression">
            <a
              onClick={() => {
                setActiveNav(2);
              }}
              className={`${styles.contentLink} ${
                activeNav === 2 || inView2 ? styles.active : ""
              }`}
            >
              Explore gene expression
            </a>
          </Link>
          <Link passHref href="#analyze-datasets">
            <a
              onClick={() => {
                setActiveNav(3);
              }}
              className={`${styles.contentLink} ${
                activeNav === 3 || inView3 ? styles.active : ""
              }`}
            >
              Analyze datasets
            </a>
          </Link>
          <Link passHref href="#download-data">
            <a
              onClick={() => {
                setActiveNav(4);
              }}
              className={`${styles.contentLink} ${
                activeNav === 4 || inView4 ? styles.active : ""
              }`}
            >
              Download data
            </a>
          </Link>
          <Link passHref href="#expedite-collaboration">
            <a
              onClick={() => {
                setActiveNav(5);
              }}
              className={`${styles.contentLink} ${
                activeNav === 5 || inView5 ? styles.active : ""
              }`}
            >
              Expedite collaboration
            </a>
          </Link>
        </div>
        <div className={styles.contentContainer}>
          <div className={styles.contentRow} id="single-cell" ref={section1}>
            <div
              className={`${styles.contentInfoCol} ${
                inView1 ? styles.active : null
              }`}
            >
              <div className={styles.contentInfoFigureCol}>
                <div className={styles.figureWrapper}>
                  <SingleCellDataIcon />
                </div>
                <span className={styles.figureSeparator}></span>
              </div>
              <div className={styles.contentInfoTextCol}>
                <h2>Quickly find the single cell data you need</h2>
                <p>
                  Browse hundreds of standardized data collections and millions
                  of cells characterizing the functionality of healthy mouse and
                  human tissues.
                </p>
                <div className={styles.linkContainer}>
                  <Link href="/" passHref>
                    <a>
                      Browse data collections
                      <span className={styles.linkArrow}>
                        <LinkArrow />
                      </span>
                    </a>
                  </Link>
                </div>
              </div>
            </div>
            <div className={styles.contentImageCol}>
              <Image src={SingleCellDataImg} alt="" />
            </div>
          </div>

          <div
            className={styles.contentRow}
            id="gene-expression"
            ref={section2}
          >
            <div
              className={`${styles.contentInfoCol} ${
                inView2 ? styles.active : ""
              }`}
            >
              <div className={styles.contentInfoFigureCol}>
                <div className={styles.figureWrapper}>
                  <GeneExpressionIcon />
                </div>
                <span className={styles.figureSeparator}></span>
              </div>
              <div className={styles.contentInfoTextCol}>
                <h2>Explore gene expression across tissues and cell types</h2>
                <p>
                  Visualize the expression of genes and gene sets using the
                  largest integrated resource of 10s of millions of cells.
                </p>
                <div className={styles.linkContainer}>
                  <Link href="/" passHref>
                    <a>
                      See how it works
                      <span className={styles.linkArrow}>
                        <LinkArrow />
                      </span>
                    </a>
                  </Link>
                  <Link href="/" passHref>
                    <a>
                      Find my gene
                      <span className={styles.linkArrow}>
                        <LinkArrow />
                      </span>
                    </a>
                  </Link>
                </div>
              </div>
            </div>
            <div className={styles.contentImageCol}>
              <Image src={GeneExpressionImg} alt="" />
            </div>
          </div>

          <div
            className={styles.contentRow}
            id="analyze-datasets"
            ref={section3}
          >
            <div
              className={`${styles.contentInfoCol} ${
                inView3 ? styles.active : ""
              }`}
            >
              <div className={styles.contentInfoFigureCol}>
                <div className={styles.figureWrapper}>
                  <AnalyzeDatasetsIcon />
                </div>
                <span className={styles.figureSeparator}></span>
              </div>
              <div className={styles.contentInfoTextCol}>
                <h2>
                  Execute on-demand interactive analyses of single cell datasets
                </h2>
                <p>
                  Visually explore how patterns of gene expression are
                  determined by environmental and genetic factors using an
                  interactive speed no-code UI. Understand published datasets or
                  use them as a launchpad to identify new cell sub-types and
                  states.
                </p>
                <div className={styles.linkContainer}>
                  <Link href="/" passHref>
                    <a>
                      Explore a multi-tissue atlas
                      <span className={styles.linkArrow}>
                        <LinkArrow />
                      </span>
                    </a>
                  </Link>
                  <Link href="/" passHref>
                    <a>
                      See all the studies you can explore
                      <span className={styles.linkArrow}>
                        <LinkArrow />
                      </span>
                    </a>
                  </Link>
                </div>
              </div>
            </div>
            <div className={styles.contentImageCol}>
              <Image src={AnalyzeDatasetsImg} alt="" />
            </div>
          </div>

          <div className={styles.contentRow} id="download-data" ref={section4}>
            <div
              className={`${styles.contentInfoCol} ${
                inView4 ? styles.active : ""
              }`}
            >
              <div className={styles.contentInfoFigureCol}>
                <div className={styles.figureWrapper}>
                  <DownloadDataIcon />
                </div>
                <span className={styles.figureSeparator}></span>
              </div>
              <div className={styles.contentInfoTextCol}>
                <h2>Download and integrate data with zero wrangling</h2>
                <p>
                  Integrate datasets with zero data wrangling. Datasets with
                  standard metadata annotations can be downloaded in AnnData and
                  Seurat formats and custom cell selections from the complete
                  corpus can be downloaded directly from R and python.
                </p>
                <div className={styles.linkContainer}>
                  <Link href="/" passHref>
                    <a>
                      Browse datasets for download
                      <span className={styles.linkArrow}>
                        <LinkArrow />
                      </span>
                    </a>
                  </Link>
                  <span>
                    Download with <a href="/">R</a> or <a href="/">Python</a>
                    <span className={styles.linkArrow}>
                      <LinkArrow />
                    </span>
                  </span>
                </div>
              </div>
            </div>
            <div className={styles.contentImageCol}>
              <Image src={DownloadDataImg} alt="" />
            </div>
          </div>

          <div
            className={styles.contentRow}
            id="expedite-collaboration"
            ref={section5}
          >
            <div
              className={`${styles.contentInfoCol} ${
                inView5 ? styles.active : ""
              }`}
            >
              <div className={styles.contentInfoFigureCol}>
                <div className={styles.figureWrapper}>
                  <ExpediteCollaborationIcon />
                </div>
              </div>
              <div className={styles.contentInfoTextCol}>
                <h2>Expedite collaborative data analysis</h2>
                <p>
                  Eliminate communication overhead and expedite cell type
                  characterization by empowering tissue experts to directly
                  explore and annotate datasets.
                </p>
                <div className={styles.linkContainer}>
                  <Link href="/" passHref>
                    <a>
                      Learn how it works
                      <span className={styles.linkArrow}>
                        <LinkArrow />
                      </span>
                    </a>
                  </Link>
                </div>
              </div>
            </div>
            <div className={styles.contentImageCol}>
              <Image src={ExpediteCollaborationImg} alt="" />
            </div>
          </div>
        </div>

        <div className={styles.publications}>
          <div className={styles.sciencePublications}>
            <span className={styles.pubSectionTitle}>Publications</span>

            <div className={styles.pubArticle}>
              <span className={styles.pubArticleDate}>03.03.21 - medRxiv</span>
              <h3 className={styles.pubArticleTitle}>
                Peng J et al.: Estimation of secondary household attack rates
                for emergent SARS-CoV-2 variants detected by genomic
                surveillance at a community-based testing site in San Francisco.
              </h3>
              <p className={styles.pubArticleCitation}>
                medRxiv, 2021 doi:10.1101/2021.03.01.21252705
              </p>
              <a className={styles.pubArticleLink} href="#">
                Read More
                <span className={styles.linkArrow}>
                  <LinkArrow />
                </span>
              </a>
            </div>

            <div className={styles.pubArticle}>
              <span className={styles.pubArticleDate}>03.03.21 - medRxiv</span>
              <h3 className={styles.pubArticleTitle}>
                Peng J et al.: Estimation of secondary household attack rates
                for emergent SARS-CoV-2 variants detected by genomic
                surveillance at a community-based testing site in San Francisco.
              </h3>
              <p className={styles.pubArticleCitation}>
                medRxiv, 2021 doi:10.1101/2021.03.01.21252705
              </p>
              <a className={styles.pubArticleLink} href="#">
                Read More
                <span className={styles.linkArrow}>
                  <LinkArrow />
                </span>
              </a>
            </div>

            <div className={styles.pubArticle}>
              <span className={styles.pubArticleDate}>03.03.21 - medRxiv</span>
              <h3 className={styles.pubArticleTitle}>
                Peng J et al.: Estimation of secondary household attack rates
                for emergent SARS-CoV-2 variants detected by genomic
                surveillance at a community-based testing site in San Francisco.
              </h3>
              <p className={styles.pubArticleCitation}>
                medRxiv, 2021 doi:10.1101/2021.03.01.21252705
              </p>
              <a className={styles.pubArticleLink} href="#">
                Read More
                <span className={styles.linkArrow}>
                  <LinkArrow />
                </span>
              </a>
            </div>
          </div>
          <div className={styles.publicationsSeparator}></div>
          <div className={styles.newsPublications}>
            <span className={styles.pubSectionTitle}>
              CELL X GENE IN THE NEWS
            </span>
            <div className={styles.pubSectionImage}>
              <Image src={NewsImage1} alt="" />
            </div>
            <div className={styles.pubSectionImage}>
              <Image src={NewsImage2} alt="" />
            </div>
          </div>
        </div>
      </main>
      <footer className={styles.footer}>
        <div className={styles.footerTopContainer}>
          <div className={styles.footerLogo}>
            <Link href="/" passHref>
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
            <a href="#">Github</a>
            <a href="#">Careers</a>
            <a href="#">Resources</a>
          </div>
        </div>

        <div className={styles.footerBottomContainer}>
          <div className={styles.footerBottomLinks}>
            <a href="#">Privacy</a>
            <a href="#">Terms</a>
            <a href="#">Contact Us</a>
          </div>

          <div className={styles.footerBottomLogos}>
            <span className={styles.footerBottomLogoText}>
              In partnership with:
            </span>
            <div className={styles.footerBottomLogosInner}>
              <a href="#" className={styles.footerBottomLogoLeft}>
                <Image src={CZILogo} alt="Chan Zuckerberg Initiative logo" />
              </a>
              <span className={styles.footerBottomLogoSeparator}></span>
              <a href="#" className={styles.footerBottomLogoRight}>
                <Image src={BiohubLogo} alt="CZ Biohub logo" />
              </a>
            </div>
          </div>
        </div>
      </footer>
    </>
  );
};

export default LandingPage;
