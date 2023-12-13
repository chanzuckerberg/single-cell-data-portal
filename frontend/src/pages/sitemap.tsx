import styled from "@emotion/styled";
import fs from "fs";
import { GetServerSideProps } from "next";
import Head from "next/head";
import pathTool from "path";
import { ROUTES } from "src/common/constants/routes";
import {
  useCellTypeMetadata,
  useTissueMetadata,
} from "src/common/queries/cellGuide";
import { useMemo, useState, useEffect } from "react";
import { fontWeightSemibold } from "src/common/theme";
import { MAX_WIDTH_BREAKPOINT_PX } from "src/components/MobileFriendlyHeader/style";

const DOC_SITE_FOLDER_NAME = "doc-site";

const DownChevron = () => (
  <svg viewBox="0 0 14 9" fill="none" xmlns="http://www.w3.org/2000/svg">
    <path
      d="M12.3136 1.00012L6.65674 6.65697L0.999884 1.00012"
      stroke="black"
      strokeWidth="2"
    />
  </svg>
);

const AccordionOpen = () => (
  <svg viewBox="0 0 42 42" fill="none" xmlns="http://www.w3.org/2000/svg">
    <circle
      cx="21"
      cy="21"
      r="21"
      transform="matrix(1 0 0 -1 0 42)"
      fill="#3867FA"
    />
    <path
      d="M26.6569 22.9999L21 17.343L15.3431 22.9999"
      stroke="white"
      strokeWidth="2"
    />
  </svg>
);

const AccordionClosed = () => (
  <svg viewBox="0 0 42 43" fill="none" xmlns="http://www.w3.org/2000/svg">
    <circle
      cx="21"
      cy="21.6799"
      r="20.5512"
      stroke="#3867FA"
      strokeWidth="0.897512"
    />
    <path
      d="M26.6569 19.68L21 25.3369L15.3431 19.68"
      stroke="#3867FA"
      strokeWidth="2"
    />
  </svg>
);

const SitemapLayout = styled.div`
  max-width: 1400px;
  margin: auto;
  padding-top: 80px;
  padding-bottom: 120px;
  padding-left: 120px;
  padding-right: 120px;

  @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
    padding: 40px 25px;
  }
`;

const SitemapTitle = styled.h1`
  margin-bottom: 0;
  font-size: 42px;
  line-height: 56.7px;
  font-weight: ${fontWeightSemibold};

  @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
    font-size: 28px;
  }
`;

const SitemapNav = styled.nav`
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(250px, 1fr));
  white-space: nowrap;

  @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
    grid-template-columns: repeat(auto-fill, 50%);
    margin-top: 10px;
  }
`;

const SitemapPageLink = styled.a`
  display: inline-block;
  margin-top: 40px;
  font-size: 22px;
  line-height: 28px;
  font-weight: ${fontWeightSemibold};
  color: #000000;

  &:hover {
    text-decoration: none;
    color: #000000;
  }

  @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
    font-size: 14px;
    line-height: 28px;
    margin-top: 10px;
  }

  svg {
    width: 12px;
    margin-left: 12px;

    @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
      width: 10px;
      margin-left: 10px;
    }
  }
`;

const SitemapNavLink = styled.a`
  display: inline-block;
  margin-top: 40px;
  font-size: 22px;
  line-height: 28px;
  font-weight: ${fontWeightSemibold};
  color: #0073ff;

  &:hover {
    text-decoration: none;
    color: #0056c6;
  }

  @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
    font-size: 14px;
    line-height: 28px;
    margin-top: 10px;
  }
`;

const SitemapSection = styled.section`
  border-top: 3px solid #f1f1f1;
  padding-top: 75px;
  margin-top: 75px;

  @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
    padding-top: 40px;
    margin-top: 40px;
  }

  h2,
  h2 a {
    font-size: 22px;
    font-weight: ${fontWeightSemibold};
    line-height: 28px;
    color: #0073ff;

    &:hover {
      text-decoration: none;
      color: #0056c6;
    }
  }
`;

const SitemapSectionGrid = styled.div`
  display: grid;
  grid-template-columns: 1fr 1fr;
  grid-gap: 0 60px;

  @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
    grid-template-columns: 1fr;
  }

  a {
    display: inline-block;
    margin-top: 20px;
    font-size: 16px;
    line-height: 26px;
    color: #0073ff;

    &:hover {
      text-decoration: none;
      color: #0056c6;
    }

    @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
      font-size: 14px;
    }
  }

  p {
    margin-top: 25px;
    font-size: 16px;
    line-height: 26px;
    font-weight: 500;
    color: #545454;
    margin-bottom: 0;

    @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
      font-size: 14px;
    }
  }

  li {
    &::marker {
      color: #0073ff;
    }

    &:hover {
      &::marker {
        color: #0056c6;
      }
    }
  }
`;

const AccordionTitle = styled.div`
  display: flex;
  align-items: center;
  justify-content: space-between;

  & h2 {
    margin: 0;
  }
`;

const AccordionContainer = styled.div`
  background: rgba(248, 249, 254, 1);
  margin-bottom: 20px;
  padding: 29px 60px 29px 86px;

  @media (max-width: 768px) {
    padding: 18px;
    margin-bottom: 15px;
  }

  &.open {
    padding: 43px 60px 50px 86px;

    @media (max-width: 768px) {
      padding: 22px 18px 40px 18px;
    }
  }
`;

const AccordionIconContainer = styled.div`
  width: 40px;
  height: 40px;
  cursor: pointer;

  @media (max-width: 768px) {
    width: 24px;
    height: 24px;
  }

  & svg {
    width: 100%;
    height: 100%;
  }
`;

interface Directory {
  dirName: string;
  files: Array<string>;
  slug: Array<string>;
  subDirectories: Array<Directory>;
}

const CACHED_FILE_PATHS = new Map<string, Directory>();
function getDocPaths(...root: Array<string>): Directory {
  const cacheStore = CACHED_FILE_PATHS && CACHED_FILE_PATHS.get(root.join("/"));
  if (cacheStore) {
    return cacheStore;
  }

  const newDirectory = {
    dirName: root[root.length - 1] ?? "",
    files: [],
    slug: root,
    subDirectories: [],
  } as Directory;

  const paths = fs.readdirSync(pathTool.join(DOC_SITE_FOLDER_NAME, ...root));
  paths.forEach((path) => {
    const builtFullPath = pathTool.join(DOC_SITE_FOLDER_NAME, ...root, path);
    if (fs.lstatSync(builtFullPath).isFile()) {
      if (path.endsWith(".mdx"))
        newDirectory.files.push(path.replace(".mdx", ""));
    } else {
      if (getDocPaths(...root, path).dirName !== ".dev") {
        newDirectory.subDirectories.push(getDocPaths(...root, path));
      }
    }
  });
  CACHED_FILE_PATHS.set(root.join("/"), newDirectory);
  return newDirectory;
}

export const getServerSideProps: GetServerSideProps = async () => {
  const response = await fetch(
    "https://api.cellxgene.cziscience.com/dp/v1/collections/index"
  );

  const collections: Collection[] = await response.json();
  const docPaths = getDocPaths();

  return {
    props: {
      collections,
      docPaths,
    },
  };
};

interface Props {
  docPaths: {
    files: [];
    subDirectories: [];
  };
  collections: Collection[];
  cellTypeMetadata: CellTypeMetadata[];
}

interface Directory {
  dirName: string;
  files: Array<string>;
  subDirectories: Directory[];
}

interface Collection {
  id: string;
  name: string;
}

interface CellTypeMetadata {
  id: string;
  name: string;
}

const Sitemap = ({ docPaths, collections }: Props): JSX.Element => {
  const { files, subDirectories } = docPaths;

  const { data: cellTypes } = useCellTypeMetadata();
  const { data: tissueData } = useTissueMetadata();

  const cellGuideData: CellTypeMetadata[] = useMemo(() => {
    if (!cellTypes || !tissueData) return [];
    const entities: CellTypeMetadata[] = [];
    for (const cellType in cellTypes) {
      entities.push({
        ...cellTypes[cellType],
      });
    }
    for (const tissue in tissueData) {
      entities.push({
        ...tissueData[tissue],
      });
    }
    return entities.sort((a, b) =>
      String(a.name).localeCompare(String(b.name), undefined, {
        sensitivity: "base",
      })
    );
  }, [cellTypes, tissueData]);

  function splitMetadataByInitialLetter(sortedArray: CellTypeMetadata[]) {
    const result = Object({});

    for (const { name, id } of sortedArray) {
      const initialLetter = String(name)[0].toUpperCase();
      if (!result[initialLetter]) {
        result[initialLetter] = [];
      }
      result[initialLetter].push({ name, id });
    }

    return result;
  }

  const cellGuideMap = splitMetadataByInitialLetter(cellGuideData);

  const cellGuideLetters = Object.keys(cellGuideMap);

  const [cgAccordionsOpen, setCgAccordionsOpen] = useState<string[]>([]);

  useEffect(() => {
    if (window.innerWidth > 768) {
      setCgAccordionsOpen([cellGuideLetters[0]]);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [cellTypes]);

  return (
    <SitemapLayout>
      <Head>
        <title>Website Sitemap - CZ CELLxGENE Discover</title>
      </Head>
      <header>
        <SitemapTitle>Sitemap</SitemapTitle>
        <SitemapNav>
          <SitemapNavLink href={ROUTES.HOMEPAGE}>Home</SitemapNavLink>
          <SitemapPageLink href={`#collections`}>
            Collections
            <DownChevron />
          </SitemapPageLink>
          <SitemapNavLink href={ROUTES.DATASETS}>Datasets</SitemapNavLink>
          <SitemapNavLink href={ROUTES.WHERE_IS_MY_GENE}>
            Gene-Expression
          </SitemapNavLink>
          <SitemapPageLink href={`#docs`}>
            Docs
            <DownChevron />
          </SitemapPageLink>
          <SitemapPageLink href={`#cellguide`}>
            CellGuide
            <DownChevron />
          </SitemapPageLink>
          <SitemapNavLink href={ROUTES.PRIVACY}>Privacy</SitemapNavLink>
          <SitemapNavLink href={ROUTES.TOS}>TOS</SitemapNavLink>
        </SitemapNav>
      </header>
      <main>
        <SitemapSection id="collections">
          <h2>
            <a href={ROUTES.COLLECTIONS}>Collections</a>
          </h2>
          <SitemapSectionGrid>
            {collections.map((collection: Collection, index) => (
              <a
                href={ROUTES.COLLECTION.replace(":id", collection.id)}
                key={`collection-${index}`}
              >
                {collection.name}
              </a>
            ))}
          </SitemapSectionGrid>
        </SitemapSection>
        <SitemapSection id="docs">
          <h2>
            <a href={ROUTES.DOCS}>Docs</a>
          </h2>
          <SitemapSectionGrid>
            {files.map((file: string, index) => (
              <a key={`file-${index}`} href={`${ROUTES.DOCS}/${file}`}>
                {file.split("__")[1]}
              </a>
            ))}
            {subDirectories.map((dir: Directory, index) => (
              <div key={`dir-${index}`}>
                <p>{dir.dirName.split("__")[1]}</p>
                {dir.files.length || dir.subDirectories.length ? (
                  <ul>
                    {dir.files.map((file: string, index) => (
                      <li key={`dirFile-${index}`}>
                        <a href={`${ROUTES.DOCS}/${file}`}>
                          {file.split("__")[1]}
                        </a>
                      </li>
                    ))}
                  </ul>
                ) : null}
                <div>
                  {dir.subDirectories.map((subDir: Directory, index) => (
                    <div key={`subDir-${index}`}>
                      <p>{subDir.dirName.split("__")[1]}</p>
                      <ul>
                        {subDir.files.map((file: string, index) => (
                          <li key={`subDirFile-${index}`}>
                            <a href={`${ROUTES.DOCS}/${file}`}>
                              {file.split("__")[1]}
                            </a>
                          </li>
                        ))}
                      </ul>
                    </div>
                  ))}
                </div>
              </div>
            ))}
          </SitemapSectionGrid>
        </SitemapSection>
        <SitemapSection id="cellguide">
          <h2>
            <a href={ROUTES.CELL_GUIDE}>CellGuide</a>
          </h2>
          {cellGuideLetters.map((letter, index) => {
            const open = cgAccordionsOpen.includes(letter);

            return (
              <AccordionContainer
                className={open ? "open" : ""}
                key={`accordionContainer-${index}`}
              >
                <AccordionTitle>
                  <h2>{letter}</h2>
                  <div
                    onClick={() => {
                      let newAccordionState;

                      if (open) {
                        newAccordionState = cgAccordionsOpen.filter(
                          (a) => a !== letter
                        );

                        setCgAccordionsOpen(newAccordionState);
                      } else {
                        newAccordionState = [...cgAccordionsOpen, letter];

                        setCgAccordionsOpen(newAccordionState);
                      }
                    }}
                  >
                    <AccordionIconContainer>
                      {open ? <AccordionOpen /> : <AccordionClosed />}
                    </AccordionIconContainer>
                  </div>
                </AccordionTitle>

                {open && (
                  <SitemapSectionGrid>
                    {cellGuideMap[letter].map(
                      (item: CellTypeMetadata, index: string) => (
                        <a
                          href={
                            item.id.includes("UBERON")
                              ? `${ROUTES.CELL_GUIDE}/tissues/${item.id}`
                              : `${ROUTES.CELL_GUIDE}/${item.id}`
                          }
                          key={`cellGuideLink-${index}`}
                        >
                          {item.name}
                        </a>
                      )
                    )}
                  </SitemapSectionGrid>
                )}
              </AccordionContainer>
            );
          })}
        </SitemapSection>
      </main>
    </SitemapLayout>
  );
};

export default Sitemap;
