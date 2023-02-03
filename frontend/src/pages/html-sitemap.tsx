import fs from "fs";
import { GetServerSideProps } from "next";
import Head from "next/head";
import pathTool from "path";
import { ROUTES } from "src/common/constants/routes";
import { CommonStyle, Layout } from "src/components/common/staticPages/style";

const DOC_SITE_FOLDER_NAME = "doc-site";

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

  const collections: any[] = await response.json();
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
  collections: [];
}

interface Directory {
  dirName: string;
  files: Array<string>;
  subDirectories: Directory[];
}

const Sitemap = ({ docPaths, collections }: Props): JSX.Element => {
  const { files, subDirectories } = docPaths;

  console.log(collections);
  return (
    <Layout>
      <CommonStyle>
        <Head>
          <title>CELL&times;GENE | Sitemap</title>
        </Head>
        <header>
          <h1>Sitemap</h1>
          <nav>
            <a href={ROUTES.HOMEPAGE}>Home</a>
            <a href={ROUTES.COLLECTIONS}>Collections</a>
            <a href={ROUTES.DATASETS}>Datasets</a>
            <a href={ROUTES.WHERE_IS_MY_GENE}>Gene-Expression</a>
            <a href={ROUTES.DOCS}>Docs</a>
            <a href={ROUTES.PRIVACY}>Privacy</a>
            <a href={ROUTES.TOS}>TOS</a>
          </nav>
        </header>

        <hr />

        <main>
          <section className="sitemapSection_collections"></section>
          <section className="sitemapSection_docs">
            <div>
              files:
              <br />
              {files.map((file: String, index) => (
                <a key={`file-${index}`} href={`/${ROUTES.DOCS}/${file}`}>
                  {file.split("__")[1]}
                </a>
              ))}
            </div>
            <br />
            <div>
              subfolders:
              <br />
              {subDirectories.map((dir: Directory, index) => (
                <div key={`dir-${index}`}>
                  <p>{dir.dirName.split("__")[1]}</p>
                  {dir.files.length || dir.subDirectories.length ? (
                    <ul>
                      {dir.files.map((file: String, index) => (
                        <li key={`dirFile-${index}`}>
                          <a href={`/${ROUTES.DOCS}/${file}`}>
                            {file.split("__")[1]}
                          </a>
                        </li>
                      ))}
                      {console.log(dir.files)}
                    </ul>
                  ) : null}
                  <div>
                    {dir.subDirectories.map((subDir: Directory, index) => (
                      <div key={`subDir-${index}`}>
                        <p>{subDir.dirName.split("__")[1]}</p>
                        <ul>
                          {subDir.files.map((file: String, index) => (
                            <li key={`subDirFile-${index}`}>
                              <a href={`/${ROUTES.DOCS}/${file}`}>
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
            </div>
          </section>
        </main>
      </CommonStyle>
    </Layout>
  );
};

export default Sitemap;
