import fs from "fs";
import matter from "gray-matter";
import { GetStaticPaths } from "next";
import { MDXRemote, MDXRemoteSerializeResult } from "next-mdx-remote";
import { serialize } from "next-mdx-remote/serialize";
import Image, { ImageProps } from "next/image";
import NextLink from "next/link";
import pathTool from "path";
import { Fragment, memo } from "react";
import EmbeddedGoogleSlides from "src/components/EmbeddedGoogleSlides";
import Layout from "src/components/Layout";
import { StyledDocsLayout } from "src/components/Layout/style";
import styled from "styled-components";

const DOC_SITE_FOLDER_NAME = "doc-site";

interface Directory {
  dirName: string;
  files: Array<string>;
  slug: Array<string>;
  subDirectories: Array<Directory>;
}

const CACHED_FILE_PATHS = new Map<string, Directory>();
function filePaths(...root: Array<string>): Directory {
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
      newDirectory.subDirectories.push(filePaths(...root, path));
    }
  });
  CACHED_FILE_PATHS.set(root.join("/"), newDirectory);
  return newDirectory;
}

function generatePaths(
  directory: Directory,
  slugs = new Array<{ params: { slug: Array<string> } }>()
): Array<{ params: { slug: Array<string> } }> {
  directory.files.forEach((file) => {
    slugs.push({ params: { slug: [...directory.slug, file] } });
  });
  if (directory.subDirectories.length > 0) {
    directory.subDirectories.forEach((subDirectory) => {
      slugs.push(...generatePaths(subDirectory, slugs));
    });
  }
  return slugs;
}

export const getStaticPaths: GetStaticPaths = async () => {
  const filePath = filePaths();
  const paths = generatePaths(filePath);

  return {
    fallback: false,
    paths,
  };
};

export const getStaticProps = async ({
  params: { slug },
}: {
  params: { slug: Array<string> };
}): Promise<{
  props: Props;
}> => {
  const markdownWithMeta = fs.readFileSync(
    pathTool.join("doc-site", slug.join("/") + ".mdx"),
    "utf-8"
  );

  const filePath = filePaths();
  const activeFile = slug[slug.length - 1];
  const { data: frontMatter, content } = matter(markdownWithMeta);
  const mdxSource = await serialize(content);
  return {
    props: {
      activeFile,
      filePath,
      frontMatter,
      mdxSource,
      slug,
    },
  };
};

const StyledUL = styled.ul<{ $isChild: boolean }>`
  border-left: ${(props) => props.$isChild && "1px solid #ccc"};
  list-style: none;
  margin-left: ${(props) => (!props.$isChild ? "40px" : "0px")};
  margin-bottom: 0px;
  color: ${(props) => props.$isChild && "#545454"};
  width: auto;
  overflow-x: hidden;
  text-overflow: ellipsis;
  & a {
    color: inherit;
    text-decoration: inherit;
  }

  & li {
    height: 36px;
    padding: 8px 16px;
    margin-bottom: 0;
  }

  & li.active-file {
    background-color: #ffffff;
    color: #e9429a;
  }

  & li:hover {
    background-color: #eaeaea;
  }
  & div {
    padding: 0px 16px;
  }
`;
const Directory = memo(function RenderDirectory({
  activeFile,
  directory,
  isChild = false,
}: {
  activeFile: string;
  directory: Directory;
  isChild?: boolean;
}) {
  const fileComponents: Array<[string, JSX.Element]> = directory.files.map(
    (file) => {
      let href = "/docs/";
      if (directory.slug.length > 0) href += directory.slug.join("/") + "/";
      href += file;
      const isActiveFile = file === activeFile;
      return [
        file,
        <li key={file} className={isActiveFile ? "active-file" : ""}>
          <NextLink href={href} passHref>
            {file}
          </NextLink>
        </li>,
      ];
    }
  );
  const directoryComponents: Array<[string, JSX.Element]> =
    directory.subDirectories.map((directory) => {
      return [
        directory.dirName,
        <Fragment key={directory.dirName}>
          <li>{directory.dirName}</li>
          <div>
            <Directory directory={directory} isChild activeFile={activeFile} />
          </div>
        </Fragment>,
      ];
    });
  return (
    <StyledUL $isChild={isChild}>
      {[...fileComponents, ...directoryComponents]
        .sort((a, b) => {
          return a[0] < b[0] ? -1 : 1;
        })
        .map((v) => v[1])}
    </StyledUL>
  );
});

interface Props {
  activeFile: string;
  frontMatter: Record<string, any>;
  mdxSource: MDXRemoteSerializeResult;
  slug: Array<string>;
  filePath: Directory;
}

const StyledLeftNav = styled.div`
  background-color: #f8f8f8;
  border-right: 1px solid #eaeaea;
  grid-area: leftsidebar;
  width: 100%;
`;

const PageNavigator = ({
  filePath,
  activeFile,
}: {
  filePath: Directory;
  activeFile: string;
}) => {
  return (
    <StyledLeftNav>
      <Directory directory={filePath} activeFile={activeFile} />
    </StyledLeftNav>
  );
};

const DocContent = styled.div`
  width: 100%;
  overflow-x: hidden;
  overflow-wrap: break-word;
  display: flex;
  align-items: center;
  flex-direction: column;
  grid-area: content;
  & > * {
    max-width: 570px;
  }
`;

const StyledImage = styled(Image)``;

const ImageContainer = styled.div`
  width: 100%;

  > div {
    position: unset !important;
  }

  ${StyledImage} {
    object-fit: contain;
    width: 100% !important;
    position: relative !important;
    height: unset !important;
  }
`;

const DocsImage = ({ src }: ImageProps) => {
  return (
    <ImageContainer>
      <StyledImage src={src} layout={"fill"} />
    </ImageContainer>
  );
};

const components = {
  EmbeddedGoogleSlides,
  Image: DocsImage,
};
const DocPage = ({ activeFile, mdxSource, filePath }: Props) => {
  return (
    <>
      <PageNavigator filePath={filePath} activeFile={activeFile} />
      <DocContent>
        <div>
          <MDXRemote {...mdxSource} components={components} />
        </div>
      </DocContent>
    </>
  );
};

DocPage.Layout = function DocLayout({
  children,
}: {
  children: JSX.Element;
}): JSX.Element {
  return (
    <Layout>
      <StyledDocsLayout>
        <main>{children}</main>
      </StyledDocsLayout>
    </Layout>
  );
};
export default DocPage;
