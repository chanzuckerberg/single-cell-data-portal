import Head from "next/head";
import Image from "next/image";
import Link from "next/link";
import { useMemo, useRef } from "react";
import { useInView } from "react-intersection-observer";
// 4513(thuang): Comment out frameSrc for now until we figure out a compliant way to embed
// import TweetEmbed from "react-tweet-embed";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { ROUTES } from "src/common/constants/routes";
import wordmark from "src/common/images/cellxgene-discover-wordmark.svg";
import AnalyzeDatasetsImg from "src/components/common/staticPages/analyze-datasets.jpg";
import LaptopImg from "src/components/common/staticPages/cellxgene-laptop-updated.png";
import HeroBg from "src/components/common/staticPages/cellxgene_hero_bg.png";
import DownloadDataImg from "src/components/common/staticPages/download-data.png";
import ExpediteCollaborationImg from "src/components/common/staticPages/expedite-collaborative-analysis.png";
import FreethinkThumbnailImg from "src/components/common/staticPages/freethink-video-thumbnail.png";
import GeneExpressionImg from "src/components/common/staticPages/explore-gene-expression.png";
import SingleCellDataImg from "src/components/common/staticPages/single-cell-data.png";
import AnalyzeDatasetsIcon from "./icons/analyze-datasets";
import CellxgeneIconSmall from "./icons/cellxgene-icon-small";
import DownloadDataIcon from "./icons/download-data";
import ExpediteCollaborationIcon from "./icons/expedite-collaboration";
import LinkArrow from "./icons/external-link-arrow";
import GeneExpressionIcon from "./icons/gene-expression";
import SingleCellDataIconActive from "./icons/single-cell-data-active";
import SingleCellDataIconInactive from "./icons/single-cell-data-inactive";
import styles from "./index.module.scss";
import { useViewMode } from "src/common/hooks/useViewMode";
import { useFetchDatasets } from "src/common/queries/filter";
import {
  LANDING_PAGE_CELLS_HERO_NUM_ID,
  LANDING_PAGE_CELLTYPES_HERO_NUM_ID,
  LANDING_PAGE_DATASETS_HERO_NUM_ID,
  LANDING_PAGE_FALLBACK_CELLS_HERO_NUM,
  LANDING_PAGE_FALLBACK_CELLTYPES_HERO_NUM,
  LANDING_PAGE_FALLBACK_DATASETS_HERO_NUM,
} from "./constants";

const ROOT_MARGIN = "-50% 0px -50% 0px";

const LandingPage = (): JSX.Element => {
  const { ref: observerSection1, inView: inView1 } = useInView({
    rootMargin: ROOT_MARGIN,
  });
  const scrollSection1 = useRef<HTMLDivElement>(null);

  const { ref: observerSection2, inView: inView2 } = useInView({
    rootMargin: ROOT_MARGIN,
  });
  const scrollSection2 = useRef<HTMLDivElement>(null);

  const { ref: observerSection3, inView: inView3 } = useInView({
    rootMargin: ROOT_MARGIN,
  });
  const scrollSection3 = useRef<HTMLDivElement>(null);

  const { ref: observerSection4, inView: inView4 } = useInView({
    rootMargin: ROOT_MARGIN,
  });
  const scrollSection4 = useRef<HTMLDivElement>(null);

  const { ref: observerSection5, inView: inView5 } = useInView({
    rootMargin: ROOT_MARGIN,
  });
  const scrollSection5 = useRef<HTMLDivElement>(null);

  const { mode, status } = useViewMode();
  const { data, isLoading, isSuccess } = useFetchDatasets(mode, status);

  const cellsHeroNum: string | null = useMemo(() => {
    if (isLoading) return null;
    if (!data || !isSuccess) return LANDING_PAGE_FALLBACK_CELLS_HERO_NUM;
    const total = data.reduce(
      (acc, curr) => acc + (curr.primary_cell_count ?? 0),
      0
    );

    if (total === 0) return LANDING_PAGE_FALLBACK_CELLS_HERO_NUM;

    const formatter = new Intl.NumberFormat("en", {
      notation: "compact",
      compactDisplay: "short",
      maximumFractionDigits: 1,
    });
    return formatter.format(total);
  }, [data, isLoading, isSuccess]);

  const cellTypesHeroNum: string | null = useMemo(() => {
    if (isLoading) return null;
    if (!data || !isSuccess) return LANDING_PAGE_FALLBACK_CELLTYPES_HERO_NUM;
    // add the cell types to the set and then get the length of the set
    const unique_cell_types = new Set();
    data.forEach((dataset) => {
      dataset.cell_type.forEach((cell_type) => {
        unique_cell_types.add(cell_type.ontology_term_id);
      });
    });
    if (unique_cell_types.size === 0) {
      return LANDING_PAGE_FALLBACK_CELLTYPES_HERO_NUM;
    }
    return `${unique_cell_types.size}`;
  }, [data, isLoading, isSuccess]);

  const datasetsHeroNum: string | null = useMemo(() => {
    if (isLoading) return null;
    if (!data || !isSuccess || Object.keys(data).length === 0) {
      return LANDING_PAGE_FALLBACK_DATASETS_HERO_NUM;
    }
    // add the cell types to the set and then get the length of the set
    return `${Object.keys(data).length}`;
  }, [data, isLoading, isSuccess]);

  const SUB_HEADING = "13.05.22 - CZ CELLxGENE";

  const CTA_TEXT_EXPLORE_DATASETS = "Explore Datasets";

  const publications = [
    {
      title:
        "The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans",
      links: [
        {
          subheading: "03.04.22 - bioRxiv",
          ctaLink: "https://doi.org/10.1101/2021.07.19.452956",
          ctaText: "Read More",
          ctaHighlight: false,
          ctaLogo: false,
        },
        {
          subheading: "13.05.22 - Science",
          ctaLink: "https://www.science.org/doi/10.1126/science.abl4896",
          ctaText: "Read More",
          ctaHighlight: false,
          ctaLogo: false,
        },
        {
          subheading: SUB_HEADING,
          ctaLink:
            "https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5",
          ctaText: CTA_TEXT_EXPLORE_DATASETS,
          ctaHighlight: true,
          ctaLogo: true,
        },
      ],
    },
    {
      title:
        "Cross-tissue immune cell analysis reveals tissue-specific features in humans",
      links: [
        {
          subheading: "13.05.22 - Science",
          ctaLink: "https://doi.org/10.1126/science.abl5197",
          ctaText: "Read More",
          ctaHighlight: false,
          ctaLogo: false,
        },
        {
          subheading: "20.07.21 - bioRxiv",
          ctaLink: "https://doi.org/10.1101/2021.04.28.441762",
          ctaText: "Read More",
          ctaHighlight: false,
          ctaLogo: false,
        },
        {
          subheading: SUB_HEADING,
          ctaLink:
            "https://cellxgene.cziscience.com/collections/62ef75e4-cbea-454e-a0ce-998ec40223d3",
          ctaText: CTA_TEXT_EXPLORE_DATASETS,
          ctaHighlight: true,
          ctaLogo: true,
        },
      ],
    },
    {
      title: "An integrated cell atlas of the human lung in health and disease",
      links: [
        {
          subheading: "11.03.21 - bioRxiv",
          ctaLink: "https://doi.org/10.1101/2022.03.10.483747",
          ctaText: "Read More",
          ctaHighlight: false,
          ctaLogo: false,
        },
        {
          subheading: SUB_HEADING,
          ctaLink:
            "https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293",
          ctaText: CTA_TEXT_EXPLORE_DATASETS,
          ctaHighlight: true,
          ctaLogo: true,
        },
      ],
    },
    {
      title:
        "Single-cell eQTL mapping identifies cell type–specific genetic control of autoimmune disease",
      links: [
        {
          subheading: "08.04.22 - Science",
          ctaLink: "https://doi.org/10.1126/science.abf3041",
          ctaText: "Read More",
          ctaHighlight: false,
          ctaLogo: false,
        },
        {
          subheading: SUB_HEADING,
          ctaLink:
            "https://cellxgene.cziscience.com/collections/dde06e0f-ab3b-46be-96a2-a8082383c4a1",
          ctaText: CTA_TEXT_EXPLORE_DATASETS,
          ctaHighlight: true,
          ctaLogo: true,
        },
      ],
    },
    {
      title:
        "Single-cell RNA-seq reveals cell type–specific molecular and genetic associations to lupus",
      links: [
        {
          subheading: "08.04.22 - Science",
          ctaLink: "https://doi.org/10.1126/science.abf1970",
          ctaText: "Read More",
          ctaHighlight: false,
          ctaLogo: false,
        },
        {
          subheading: SUB_HEADING,
          ctaLink:
            "https://cellxgene.cziscience.com/collections/436154da-bcf1-4130-9c8b-120ff9a888f2",
          ctaText: CTA_TEXT_EXPLORE_DATASETS,
          ctaHighlight: true,
          ctaLogo: true,
        },
      ],
    },
  ];

  const newsLinks = [
    {
      title: `Accelerating COVID-19 Research with New Single-Cell Technologies`,
      link: `https://cziscience.medium.com/accelerating-covid-19-research-with-new-single-cell-technologies-7cefdfc54cbb`,
      newTab: true,
    },
    {
      title: `5 Ways Single-Cell Biology Is Advancing Our Understanding of Disease`,
      link: `https://chanzuckerberg.com/blog/single-cell-biology-understanding-disease/`,
      newTab: true,
    },
  ];

  return (
    <>
      <Head>
        <title>CZ CELLxGENE Discover - Cellular Visualization Tool</title>
        <meta
          name="description"
          content="Chan Zuckerberg CELLxGENE Discover is a tool to find, download, and visually explore curated and standardized single-cell biology datasets."
        />
        <meta
          property="og:title"
          content="Chan Zuckerberg CELLxGENE Discover"
        />
        <meta
          property="og:description"
          content="Find, download and visually explore curated and standardized single-cell datasets"
        />
        <meta property="twitter:card" content="summary" />
        <meta
          property="twitter:title"
          content="Chan Zuckerberg CELLxGENE Discover"
        />
        <meta
          property="twitter:description"
          content="Find, download and visually explore curated and standardized single-cell datasets"
        />
        <meta
          property="og:image"
          content="https://cellxgene.cziscience.com/open-graph.jpg"
        />
        <meta
          property="twitter:image"
          content="https://cellxgene.cziscience.com/open-graph.jpg"
        />
        <meta property="og:type" content="website" />
      </Head>
      <div className={styles.landingPage}>
        <div
          className={styles.heroContainer}
          style={{
            backgroundImage: `url(${HeroBg.src})`,
            backgroundPosition: "left bottom",
            backgroundRepeat: "no-repeat",
          }}
        >
          <div className={styles.heroImgContainer}>
            <div
              className={styles.laptopImg}
              data-testid="laptop-with-cell-data-on-screen"
            >
              <Image
                src={LaptopImg}
                alt="laptop with cell data on screen"
                priority
              />
            </div>
          </div>
          <div className={styles.heroTextContainer}>
            <div className={styles.heroLogo}>
              <Image src={wordmark} alt="CZ CELLxGENE Discover" priority />
            </div>
            <h1>Discover the mechanisms of human health</h1>
            <p>
              Download and visually explore reference-quality data to understand
              the functionality of human tissues at the cellular level with Chan
              Zuckerberg CELL by GENE Discover (CZ CELLxGENE Discover).
            </p>
            <div className={styles.heroStatsContainer}>
              <div>
                <span>unique cells</span>
                <p data-testid={LANDING_PAGE_CELLS_HERO_NUM_ID}>
                  {cellsHeroNum}
                </p>
              </div>
              <div>
                <span>datasets</span>
                <p data-testid={LANDING_PAGE_DATASETS_HERO_NUM_ID}>
                  {datasetsHeroNum}
                </p>
              </div>
              <div>
                <span>cell types</span>
                <p data-testid={LANDING_PAGE_CELLTYPES_HERO_NUM_ID}>
                  {cellTypesHeroNum}
                </p>
              </div>
            </div>
          </div>
        </div>
        <div className={styles.main}>
          <div>
            <div className={styles.contentNav}>
              <div className={styles.contentNavSubrow}>
                <div
                  className={`${styles.contentLink} ${
                    inView1 ? styles.active : ""
                  }`}
                  onClick={() => {
                    scrollSection1.current?.scrollIntoView({
                      behavior: "smooth",
                    });

                    track(EVENTS.HOMEPAGE_LEARN_FIND_SINGLE_CELL_DATA_CLICKED);
                  }}
                >
                  Find single-cell data
                </div>
                <div
                  className={`${styles.contentLink} ${
                    inView2 ? styles.active : ""
                  }`}
                  onClick={() => {
                    scrollSection2.current?.scrollIntoView({
                      behavior: "smooth",
                    });

                    track(
                      EVENTS.HOMEPAGE_LEARN_EXPLORE_GENE_EXPRESSION_CLICKED
                    );
                  }}
                >
                  Explore gene expression
                </div>
              </div>
              <div className={styles.contentNavSubrow}>
                <div
                  className={`${styles.contentLink} ${
                    inView3 ? styles.active : ""
                  }`}
                  onClick={() => {
                    scrollSection3.current?.scrollIntoView({
                      behavior: "smooth",
                    });

                    track(EVENTS.HOMEPAGE_LEARN_ANALYZE_DATASETS_CLICKED);
                  }}
                >
                  Analyze datasets
                </div>
                <div
                  className={`${styles.contentLink} ${
                    inView4 ? styles.active : ""
                  }`}
                  onClick={() => {
                    scrollSection4.current?.scrollIntoView({
                      behavior: "smooth",
                    });

                    track(EVENTS.HOMEPAGE_LEARN_DOWNLOAD_DATA_CLICKED);
                  }}
                >
                  Download data
                </div>
                <div
                  className={`${styles.contentLink} ${
                    inView5 ? styles.active : ""
                  }`}
                  onClick={() => {
                    scrollSection5.current?.scrollIntoView({
                      behavior: "smooth",
                    });

                    track(EVENTS.HOMEPAGE_LEARN_EXPEDITE_COLLABORATION_CLICKED);
                  }}
                >
                  Expedite collaboration
                </div>
              </div>
            </div>
            <div className={styles.contentOverflowWrapper}>
              <div className={styles.contentContainer}>
                <div ref={observerSection1}>
                  <div
                    className={`${styles.contentRow} ${styles.contentRowFirst}`}
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
                          {inView1 ? (
                            <SingleCellDataIconActive />
                          ) : (
                            <SingleCellDataIconInactive />
                          )}
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
                            <a
                              onClick={() =>
                                track(EVENTS.BROWSE_COLLECTIONS_CLICKED, {
                                  button: "browse data collections",
                                })
                              }
                            >
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
                          alt="CELLxGENE Discover collections page with sidebar for filtering through the table of data"
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
                          Visualize the expression of genes and gene sets using
                          the largest integrated resource of over 25 million
                          cells.
                        </p>
                        <div className={styles.linkContainer}>
                          <Link href={ROUTES.WHERE_IS_MY_GENE} passHref>
                            <a onClick={() => track(EVENTS.WMG_CLICKED, {})}>
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
                          determined by environmental and genetic factors using
                          an interactive speed no-code UI. Understand published
                          datasets or use them as a launchpad to identify new
                          cell sub-types and states.
                        </p>
                        <div className={styles.linkContainer}>
                          <a
                            /**
                             * (thuang): Open in a new tab, so we don't lose
                             * analytics API call that would otherwise get
                             * cancelled due to a new full page request
                             */
                            href={`${ROUTES.HOMEPAGE}e/53d208b0-2cfd-4366-9866-c3c6114081bc.cxg/`}
                            rel="noopener"
                            target="_blank"
                            onClick={() =>
                              track(EVENTS.DATASET_EXPLORE_CLICKED, {
                                // (thuang): Please update the dataset name when the href link changes
                                // to a different dataset
                                dataset_name: "Tabula Sapiens - All Cells",
                              })
                            }
                          >
                            Explore a multi-tissue atlas
                            <span className={styles.linkArrow}>
                              <LinkArrow />
                            </span>
                          </a>
                          <Link href={ROUTES.COLLECTIONS} passHref>
                            <a
                              onClick={() =>
                                track(EVENTS.BROWSE_COLLECTIONS_CLICKED, {
                                  button: "explore the studies",
                                })
                              }
                            >
                              Explore the studies
                              <span className={styles.linkArrow}>
                                <LinkArrow />
                              </span>
                            </a>
                          </Link>
                          {/* DOC PAGE LINKS NEED TO BE OPENED IN A NEW TAB IN ORDER TO LOAD UNIQUE CSP DIRECTIVE */}
                          <a
                            href={`${ROUTES.DOCS}/04__Analyze%20Public%20Data/4_1__Hosted%20Tutorials`}
                            rel="noopener"
                            target="_blank"
                            onClick={() =>
                              track(EVENTS.BROWSE_TUTORIALS_CLICKED)
                            }
                          >
                            Browse tutorials
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
                          Integrate datasets with zero data wrangling. Datasets
                          with standard metadata annotations can be downloaded
                          in AnnData and Seurat formats.
                        </p>
                        <div className={styles.linkContainer}>
                          <Link href={ROUTES.COLLECTIONS} passHref>
                            <a
                              onClick={() =>
                                track(EVENTS.BROWSE_COLLECTIONS_CLICKED, {
                                  button: "browse datasets for download",
                                })
                              }
                            >
                              Browse datasets for download
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
                          Eliminate communication overhead and expedite cell
                          type characterization by empowering tissue experts to
                          directly explore and annotate datasets.
                        </p>
                        <div className={styles.linkContainer}>
                          <a
                            onClick={() => {
                              track(
                                EVENTS.EXPLORE_CZ_CELLXGENE_ANNOTATE_CLICKED
                              );
                            }}
                            href={`${ROUTES.DOCS}/05__Annotate%20and%20Analyze%20Your%20Data/5_0__Get%20Started`}
                          >
                            Explore CZ CELLxGENE Annotate
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
          </div>

          <div className={styles.freethinkVideoSection}>
            <div className={styles.freethinkVideoSectionWrapper}>
              <div className={styles.freethinkVideoSectionTextContainer}>
                <h2>Behind-the-scenes of Chan Zuckerberg CELLxGENE</h2>
                <p>
                  Watch a behind-the-scenes look at Chan Zuckerberg CELLxGENE
                  (CZ CELLxGENE) and explore how the platform can help
                  scientists quickly surface important information that could
                  lead to discoveries in treating disease.
                </p>
              </div>
              <a
                href="https://chanzuckerberg.com/blog/behind-the-scenes-of-chan-zuckerberg-cellxgene/"
                target="_blank"
                rel="noopener"
              >
                <Image
                  src={FreethinkThumbnailImg}
                  alt="Chan Zuckerberg Initiative presents Mapping the Road to Health, a Freethink original video (external link)"
                />
              </a>
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
                      className={`${styles.pubArticleSubRow} ${
                        link.ctaHighlight
                          ? styles.pubArticleSubRowHighlight
                          : ""
                      }`}
                      key={`article-${articleIndex}-link-${linkIndex}`}
                    >
                      <span className={styles.pubArticleDate}>
                        {link.subheading}
                      </span>
                      <div className={styles.pubArticleSubRowInner}>
                        <a
                          className={styles.pubArticleLink}
                          href={link.ctaLink}
                          target="_blank"
                          rel="noopener"
                          onClick={
                            // (thuang): We use `ctaLogo` to infer a collection page link
                            // If this no longer true, we need to update this logic
                            link.ctaLogo
                              ? () =>
                                  track(EVENTS.VIEW_COLLECTION_PAGE_CLICKED, {
                                    collection_name: pub.title,
                                  })
                              : undefined
                          }
                        >
                          {link.ctaText}
                          <span className={styles.linkArrow}>
                            <LinkArrow />
                          </span>
                        </a>
                        {link.ctaLogo && (
                          <CellxgeneIconSmall
                            className={styles.pubArticleLogoIcon}
                          />
                        )}
                      </div>
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
              {/* 4513(thuang): Comment out frameSrc for now until we figure out a compliant way to embed */}
              {/* <div className={styles.pubSectionImage}>
                <TweetEmbed tweetId="1404822000464433158" />
              </div>
              <div className={styles.pubSectionImage}>
                <TweetEmbed tweetId="1396854799908282376" />
              </div> */}
              {newsLinks && (
                <div className={styles.newsLinkContainer}>
                  {newsLinks.map((link, index) => (
                    <div
                      className={styles.newsLinkItem}
                      key={`news-link-${index}`}
                    >
                      <p className={styles.newsLinkTitle}>{link.title}</p>
                      <a
                        className={styles.newsLink}
                        href={link.link}
                        target={link.newTab ? "_blank" : "_self"}
                        rel={link.newTab ? "noopener noreferrer" : ""}
                      >
                        Read More
                        <span className={styles.newsLinkArrow}>
                          <LinkArrow />
                        </span>
                      </a>
                    </div>
                  ))}
                </div>
              )}
            </div>
          </div>
        </div>
      </div>
    </>
  );
};

export default LandingPage;
