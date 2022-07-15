import { Icon } from "czifui";
import fs from "fs";
import matter from "gray-matter";
import { GetStaticPaths } from "next";
import { MDXRemote, MDXRemoteSerializeResult } from "next-mdx-remote";
import { serialize } from "next-mdx-remote/serialize";
import Image, { ImageProps } from "next/image";
import NextLink from "next/link";
import pathTool from "path";
import { Fragment, memo, useState } from "react";
import { noop } from "src/common/constants/utils";
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

const StyledUL = styled.ul<{ $isChild: boolean; isExpanded: boolean }>`
  border-left: ${(props) => props.$isChild && "1px solid #ccc"};
  list-style: none;
  margin-left: ${(props) => (!props.$isChild ? "40px" : "0px")};
  margin-bottom: 0px;
  color: ${(props) => props.$isChild && "#545454"};
  width: auto;
  overflow-x: hidden;
  text-overflow: ellipsis;
  height: ${(props) => (props.isExpanded ? "auto" : "0px")};
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

const FileListItem = ({
  file,
  isActiveFile,
  href,
}: {
  file: string;
  isActiveFile: boolean;
  href: string;
}) => {
  const formattedFileName = file.split("__")[1];

  return (
    <li key={file} className={isActiveFile ? "active-file" : ""}>
      <NextLink href={href} passHref>
        {formattedFileName}
      </NextLink>
    </li>
  );
};

const DirectoryListItem = ({
  directory,
  activeFile,
}: {
  directory: Directory;
  activeFile: string;
}) => {
  // 0 = default collapse, 1 = default expand, 2 = user collapse, 3 = user expand
  const [isExpanded, setIsExpanded] = useState<0 | 1 | 2 | 3>(0);
  return (
    <Fragment>
      <li onClick={() => setIsExpanded(isExpanded % 2 == 0 ? 3 : 2)}>
        {directory.dirName.split("__")[1]}{" "}
        <Icon
          sdsIcon={isExpanded % 2 == 1 ? "chevronDown" : "chevronRight"}
          sdsSize="s"
          sdsType="interactive"
        />
      </li>
      <div>
        <Directory
          directory={directory}
          isChild
          activeFile={activeFile}
          isExpanded={isExpanded}
          setIsExpanded={setIsExpanded}
        />
      </div>
    </Fragment>
  );
};

const Directory = memo(function RenderDirectory({
  activeFile,
  directory,
  isExpanded,
  setIsExpanded,
  isChild = false,
}: {
  activeFile: string;
  directory: Directory;
  isExpanded: 0 | 1 | 2 | 3;
  setIsExpanded: (isExpanded: 0 | 1 | 2 | 3) => void;
  isChild?: boolean;
}) {
  const fileComponents: Array<[string, JSX.Element]> = directory.files.map(
    (file) => {
      let href = "/docs/";
      if (directory.slug.length > 0) href += directory.slug.join("/") + "/";
      href += file;
      const isActiveFile = file === activeFile;
      if (isActiveFile && isExpanded < 2) setIsExpanded(1);

      return [
        file,
        <FileListItem
          key={file}
          file={file}
          href={href}
          isActiveFile={isActiveFile}
        />,
      ];
    }
  );

  const directoryComponents: Array<[string, JSX.Element | null]> =
    directory.subDirectories.map((directory) => {
      if (directory.dirName.startsWith(".")) return [directory.dirName, null];
      return [
        directory.dirName,
        <DirectoryListItem
          key={directory.dirName}
          directory={directory}
          activeFile={activeFile}
        />,
      ];
    });
  return (
    <StyledUL $isChild={isChild} isExpanded={isExpanded % 2 == 1}>
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
      <Directory
        directory={filePath}
        activeFile={activeFile}
        isExpanded={1}
        setIsExpanded={noop}
      />
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

const MDX_AVAILABLE_COMPONENTS = {
  EmbeddedGoogleSlides,
  Image: DocsImage,
};
const DocPage = ({ activeFile, mdxSource, filePath }: Props) => {
  return (
    <>
      <PageNavigator filePath={filePath} activeFile={activeFile} />
      <DocContent>
        <div>
          <MDXRemote {...mdxSource} components={MDX_AVAILABLE_COMPONENTS} />
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
