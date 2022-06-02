import Head from "next/head";
import Image from "next/image";
import Link from "next/link";
import { useRef, useState } from "react";
import { useInView } from "react-intersection-observer";
import TweetEmbed from "react-tweet-embed";
import { ROUTES } from "src/common/constants/routes";
import AnalyzeDatasetsImg from "src/components/common/staticPages/analyze-datasets.png";
import HeroBg from "src/components/common/staticPages/cellxgene_hero_bg.png";
import LaptopImg from "src/components/common/staticPages/cellxgene_laptop.png";
import DownloadDataImg from "src/components/common/staticPages/download-data.png";
import ExpediteCollaborationImg from "src/components/common/staticPages/expedite-collaboration.png";
import GeneExpressionImg from "src/components/common/staticPages/gene-expression.png";
import SingleCellDataImg from "src/components/common/staticPages/single-cell-data.png";
import AnalyzeDatasetsIcon from "./icons/analyze-datasets";
import DownloadDataIcon from "./icons/download-data";
import ExpediteCollaborationIcon from "./icons/expedite-collaboration";
import LinkArrow from "./icons/external-link-arrow";
import GeneExpressionIcon from "./icons/gene-expression";
import SingleCellDataIcon from "./icons/single-cell-data";
import styles from "./index.module.scss";

const LandingPage = (): JSX.Element => {
  const { ref: observerSection1, inView: inView1 } = useInView({
    rootMargin: "-50% 0px -50% 0px",
  });
  const scrollSection1 = useRef<HTMLDivElement>(null!);

  const { ref: observerSection2, inView: inView2 } = useInView({
    rootMargin: "-50% 0px -50% 0px",
  });
  const scrollSection2 = useRef<HTMLDivElement>(null!);

  const { ref: observerSection3, inView: inView3 } = useInView({
    rootMargin: "-50% 0px -50% 0px",
  });
  const scrollSection3 = useRef<HTMLDivElement>(null!);

  const { ref: observerSection4, inView: inView4 } = useInView({
    rootMargin: "-50% 0px -50% 0px",
  });
  const scrollSection4 = useRef<HTMLDivElement>(null!);

  const { ref: observerSection5, inView: inView5 } = useInView({
    rootMargin: "-50% 0px -50% 0px",
  });
  const scrollSection5 = useRef<HTMLDivElement>(null!);

  // HERO NUMBERS. DUMMY DATA TO BE REPLACED.
  const [cellsHeroNum] = useState("30M+");
  const [datasetsHeroNum] = useState("436");
  const [donorsHeroNum] = useState("2.7k+");

  const publications = [
    {
      title:
        "The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans",
      links: [
        {
          subheading: "03.04.22 - bioRxiv",
          ctaText: "DOI: 10.1101/2021.07.19.452956",
          ctaLink: "https://doi.org/10.1101/2021.07.19.452956",
        },
        {
          subheading: "13.05.22 - Science",
          ctaText: "DOI: 10.1126/science.abl4896",
          ctaLink: "https://www.science.org/doi/10.1126/science.abl4896",
        },
        {
          subheading: "13.05.22 - Cellxgene",
          ctaText: "The Tabula Sapiens Consortium et al. (2021) bioRxiv",
          ctaLink:
            "https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5",
        },
      ],
    },
    {
      title:
        "Cross-tissue immune cell analysis reveals tissue-specific features in humans",
      links: [
        {
          subheading: "13.05.22 - Science",
          ctaText: "DOI: 10.1126/science.abl5197",
          ctaLink: "https://doi.org/10.1126/science.abl5197",
        },
        {
          subheading: "20.07.21 - bioRxiv",
          ctaText: "DOI: 10.1101/2021.04.28.441762",
          ctaLink: "https://doi.org/10.1101/2021.04.28.441762",
        },
        {
          subheading: "13.05.22 - Cellxgene",
          ctaText: "Domínguez Conde et al. (2022) Science",
          ctaLink:
            "https://cellxgene.cziscience.com/collections/62ef75e4-cbea-454e-a0ce-998ec40223d3",
        },
      ],
    },
    {
      title: "An integrated cell atlas of the human lung in health and disease",
      links: [
        {
          subheading: "11.03.21 - bioRxiv",
          ctaText: "DOI: 10.1101/2022.03.10.483747",
          ctaLink: "https://doi.org/10.1101/2022.03.10.483747",
        },
        {
          subheading: "13.05.22 - Cellxgene",
          ctaText: "Sikkema et al. (2022) bioRxiv",
          ctaLink:
            "https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293",
        },
      ],
    },
    {
      title:
        "Single-cell eQTL mapping identifies cell type–specific genetic control of autoimmune disease",
      links: [
        {
          subheading: "08.04.22 - Science",
          ctaText: "DOI: 10.1126/science.abf3041",
          ctaLink: "https://doi.org/10.1126/science.abf3041",
        },
        {
          subheading: "13.05.22 - Cellxgene",
          ctaText: "Yazar et al. (2022) Science",
          ctaLink:
            "https://cellxgene.cziscience.com/collections/dde06e0f-ab3b-46be-96a2-a8082383c4a1",
        },
      ],
    },
    {
      title:
        "Single-cell RNA-seq reveals cell type–specific molecular and genetic associations to lupus",
      links: [
        {
          subheading: "08.04.22 - Science",
          ctaText: "DOI: 10.1126/science.abf1970",
          ctaLink: "https://doi.org/10.1126/science.abf1970",
        },
        {
          subheading: "13.05.22 - Cellxgene",
          ctaText: "Perez et al. (2022) Science",
          ctaLink:
            "https://cellxgene.cziscience.com/collections/436154da-bcf1-4130-9c8b-120ff9a888f2",
        },
      ],
    },
  ];

  return (
    <>
      <Head>
        <title>cellxgene | Home</title>
      </Head>
      <div
        className={styles.heroContainer}
        style={{
          backgroundImage: `url(${HeroBg.src})`,
          backgroundPosition: "left bottom",
          backgroundRepeat: "no-repeat",
        }}
      >
        <div className={styles.heroImgContainer}>
          <div className={styles.laptopImg}>
            <Image src={LaptopImg} alt="laptop with cell data on screen" />
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
              <p>{cellsHeroNum}</p>
            </div>
            <div>
              <span>datasets</span>
              <p>{datasetsHeroNum}</p>
            </div>
            <div>
              <span>donors</span>
              <p>{donorsHeroNum}</p>
            </div>
          </div>
        </div>
      </div>
      <div className={styles.main}>
        <div>
          <div className={styles.contentNav}>
            <div
              className={`${styles.contentLink} ${
                inView1 ? styles.active : ""
              }`}
              onClick={() => {
                scrollSection1.current?.scrollIntoView({ behavior: "smooth" });
              }}
            >
              Find single-cell data
            </div>
            <div
              className={`${styles.contentLink} ${
                inView2 ? styles.active : ""
              }`}
              onClick={() => {
                scrollSection2.current?.scrollIntoView({ behavior: "smooth" });
              }}
            >
              Explore gene expression
            </div>
            <div
              className={`${styles.contentLink} ${
                inView3 ? styles.active : ""
              }`}
              onClick={() => {
                scrollSection3.current?.scrollIntoView({ behavior: "smooth" });
              }}
            >
              Analyze datasets
            </div>
            <div
              className={`${styles.contentLink} ${
                inView4 ? styles.active : ""
              }`}
              onClick={() => {
                scrollSection4.current?.scrollIntoView({ behavior: "smooth" });
              }}
            >
              Download data
            </div>
            <div
              className={`${styles.contentLink} ${
                inView5 ? styles.active : ""
              }`}
              onClick={() => {
                scrollSection5.current?.scrollIntoView({ behavior: "smooth" });
              }}
            >
              Expedite collaboration
            </div>
          </div>
          <div className={styles.contentContainer}>
            <div ref={observerSection1}>
              <div
                className={styles.contentRow}
                id="single-cell"
                ref={scrollSection1}
              >
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
                    <h2 className={styles.mt16}>
                      Quickly find the single-cell data you need
                    </h2>
                    <p>
                      Browse hundreds of standardized data collections and
                      millions of cells characterizing the functionality of
                      healthy mouse and human tissues.
                    </p>
                    <div className={styles.linkContainer}>
                      <Link href={ROUTES.COLLECTIONS} passHref>
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
                  <div className={styles.contentImage}>
                    <Image
                      src={SingleCellDataImg}
                      alt="cellxgene collections page with sidebar for filtering through the table of data"
                    />
                  </div>
                </div>
              </div>
            </div>

            <div ref={observerSection2}>
              <div
                className={styles.contentRow}
                id="gene-expression"
                ref={scrollSection2}
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
                    <h2>
                      Explore gene expression across tissues and cell types
                    </h2>
                    <p>
                      Visualize the expression of genes and gene sets using the
                      largest integrated resource of over 30 million cells.
                    </p>
                    <div className={styles.linkContainer}>
                      <Link href={ROUTES.WHERE_IS_MY_GENE} passHref>
                        <a>
                          See how it works
                          <span className={styles.linkArrow}>
                            <LinkArrow />
                          </span>
                        </a>
                      </Link>
                    </div>
                  </div>
                </div>
                <div className={styles.contentImageCol}>
                  <div className={styles.contentImage}>
                    <Image
                      src={GeneExpressionImg}
                      alt="chart of gene expressions with genes plotted against types of tissues"
                    />
                  </div>
                </div>
              </div>
            </div>

            <div ref={observerSection3}>
              <div
                className={styles.contentRow}
                id="analyze-datasets"
                ref={scrollSection3}
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
                      Execute on-demand interactive analyses of single-cell
                      datasets
                    </h2>
                    <p>
                      Visually explore how patterns of gene expression are
                      determined by environmental and genetic factors using an
                      interactive speed no-code UI. Understand published
                      datasets or use them as a launchpad to identify new cell
                      sub-types and states.
                    </p>
                    <div className={styles.linkContainer}>
                      <Link
                        href={`${ROUTES.HOMEPAGE}/e/53d208b0-2cfd-4366-9866-c3c6114081bc.cxg/`}
                        passHref
                      >
                        <a>
                          Explore a multi-tissue atlas
                          <span className={styles.linkArrow}>
                            <LinkArrow />
                          </span>
                        </a>
                      </Link>
                      <Link href={ROUTES.COLLECTIONS} passHref>
                        <a>
                          Explore the studies
                          <span className={styles.linkArrow}>
                            <LinkArrow />
                          </span>
                        </a>
                      </Link>
                      <Link
                        href={`${ROUTES.HOMEPAGE}/docs/04_Analyze-Public-Data/4_1-Hosted-Tutorials`}
                        passHref
                      >
                        <a>
                          Browse tutorials
                          <span className={styles.linkArrow}>
                            <LinkArrow />
                          </span>
                        </a>
                      </Link>
                    </div>
                  </div>
                </div>
                <div className={styles.contentImageCol}>
                  <div className={styles.contentImage}>
                    <Image
                      src={AnalyzeDatasetsImg}
                      alt="multi-tissue visualisation with a legend showing which colors correspond to specific cell types"
                    />
                  </div>
                </div>
              </div>
            </div>

            <div ref={observerSection4}>
              <div
                className={styles.contentRow}
                id="download-data"
                ref={scrollSection4}
              >
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
                    <h2 className={styles.mt16}>
                      Download and integrate data with zero wrangling
                    </h2>
                    <p>
                      Integrate datasets with zero data wrangling. Datasets with
                      standard metadata annotations can be downloaded in AnnData
                      and Seurat formats and custom cell selections from the
                      complete corpus can be downloaded directly from R and
                      Python.
                    </p>
                    <div className={styles.linkContainer}>
                      <Link href={ROUTES.COLLECTIONS} passHref>
                        <a>
                          Browse datasets for download
                          <span className={styles.linkArrow}>
                            <LinkArrow />
                          </span>
                        </a>
                      </Link>
                      {/* LINK TO BE UPDATED POST-LAUNCH */}
                      <a href="#">
                        Download with R
                        <span className={styles.linkArrow}>
                          <LinkArrow />
                        </span>
                      </a>
                      {/* LINK TO BE UPDATED POST-LAUNCH */}
                      <a href="#">
                        Download with Python
                        <span className={styles.linkArrow}>
                          <LinkArrow />
                        </span>
                      </a>
                    </div>
                  </div>
                </div>
                <div className={styles.contentImageCol}>
                  <div className={styles.contentImage}>
                    <Image
                      src={DownloadDataImg}
                      alt="pop-up modal for downloading a dataset"
                    />
                  </div>
                </div>
              </div>
            </div>

            <div ref={observerSection5}>
              <div
                className={styles.contentRow}
                id="expedite-collaboration"
                ref={scrollSection5}
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
                    <h2 className={styles.mt16}>
                      Expedite collaborative data analysis
                    </h2>
                    <p>
                      Eliminate communication overhead and expedite cell type
                      characterization by empowering tissue experts to directly
                      explore and annotate datasets.
                    </p>
                    <div className={styles.linkContainer}>
                      {/* LINK TO BE UPDATED POST-LAUNCH */}
                      <a href="https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/desktop/quick-start.md">
                        Learn how it works
                        <span className={styles.linkArrow}>
                          <LinkArrow />
                        </span>
                      </a>
                    </div>
                  </div>
                </div>
                <div className={styles.contentImageCol}>
                  <div className={styles.contentImage}>
                    <Image
                      src={ExpediteCollaborationImg}
                      alt="pop-up modal for user to create a data directory for storing gene sets and annotations"
                    />
                  </div>
                </div>
              </div>
            </div>
          </div>
        </div>

        <div className={styles.publications}>
          <div className={styles.sciencePublications}>
            <span className={styles.pubSectionTitle}>Publications</span>

            {publications.map((pub, articleIndex) => (
              <div
                className={styles.pubArticle}
                key={`pubArticle-${articleIndex}`}
              >
                <h3 className={styles.pubArticleTitle}>{pub.title}</h3>

                {pub.links.map((link, linkIndex) => (
                  <div
                    className={styles.pubArticleSubRow}
                    key={`article-${articleIndex}-link-${linkIndex}`}
                  >
                    <div>
                      <span className={styles.pubArticleDate}>
                        {link.subheading}
                      </span>
                      <p className={styles.pubArticleCitation}>
                        {link.ctaText}
                      </p>
                    </div>
                    <a
                      className={styles.pubArticleLink}
                      href={link.ctaLink}
                      target="_blank"
                      rel="noopener"
                    >
                      Read More
                      <span className={styles.linkArrow}>
                        <LinkArrow />
                      </span>
                    </a>
                  </div>
                ))}
              </div>
            ))}
          </div>
          <div className={styles.publicationsSeparator}></div>
          <div className={styles.newsPublications}>
            <span className={`${styles.pubSectionTitle} ${styles.newsTitle}`}>
              CELL X GENE IN THE NEWS
            </span>
            <div className={styles.pubSectionImage}>
              <TweetEmbed tweetId="1404822000464433158" />
            </div>
            <div className={styles.pubSectionImage}>
              <TweetEmbed tweetId="1396854799908282376" />
            </div>
          </div>
        </div>
      </div>
    </>
  );
};

export default LandingPage;
