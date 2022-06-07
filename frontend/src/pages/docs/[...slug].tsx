import fs from "fs";
import matter from "gray-matter";
import { GetStaticPaths } from "next";
import { MDXRemote, MDXRemoteSerializeResult } from "next-mdx-remote";
import { serialize } from "next-mdx-remote/serialize";
import Image from "next/image";
import NextLink from "next/link";
import pathTool from "path";
import { Fragment, memo } from "react";
import EmbeddedGoogleSlides from "src/components/EmbeddedGoogleSlides";
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
  const { data: frontMatter, content } = matter(markdownWithMeta);
  const mdxSource = await serialize(content);
  return {
    props: {
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

  & li:hover {
    background-color: #eaeaea;
  }
  & div {
    padding: 0px 16px;
  }
`;
const Directory = memo(function RenderDirectory({
  directory,
  isChild = false,
}: {
  directory: Directory;
  isChild?: boolean;
}) {
  return (
    <StyledUL $isChild={isChild}>
      {directory.files.map((file) => {
        let href = "/docs/";
        if (directory.slug.length > 0) href += directory.slug.join("/") + "/";
        href += file;
        return (
          <li key={file}>
            <NextLink href={href} passHref>
              {file}
            </NextLink>
          </li>
        );
      })}
      {directory.subDirectories.map((directory) => {
        return (
          <Fragment key={directory.dirName}>
            <li>{directory.dirName}</li>
            <div>
              <Directory directory={directory} isChild />
            </div>
          </Fragment>
        );
      })}
    </StyledUL>
  );
});

interface Props {
  frontMatter: Record<string, any>;
  mdxSource: MDXRemoteSerializeResult;
  slug: Array<string>;
  filePath: Directory;
}

const StyledLeftNav = styled.div`
  background-color: #f8f8f8;
`;

const PageNavigator = ({ filePath }: { filePath: Directory }) => {
  return (
    <StyledLeftNav>
      <Directory directory={filePath} />
    </StyledLeftNav>
  );
};

const StyledDocsLayout = styled.div`
  display: grid;
  grid-template-columns: [nav] 368px [content] auto [tableOfContents] 368px;
`;
const DocContent = styled.div`
  width: 100%;
  overflow-x: hidden;
  overflow-wrap: break-word;
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: column;

  & > * {
    max-width: 570px;
  }
`;

const components = {
  EmbeddedGoogleSlides,
  Image,
};
const DocPage = ({ mdxSource, filePath }: Props) => {
  return (
    <StyledDocsLayout>
      <PageNavigator filePath={filePath} />
      <DocContent>
        <div>
          <MDXRemote {...mdxSource} components={components} />
        </div>
      </DocContent>
    </StyledDocsLayout>
  );
};

export default DocPage;
