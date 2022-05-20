import { Link } from "czifui";
import fs from "fs";
import matter from "gray-matter";
import { GetStaticPaths } from "next";
import { MDXRemote, MDXRemoteSerializeResult } from "next-mdx-remote";
import { serialize } from "next-mdx-remote/serialize";
import Image from "next/image";
import NextLink from "next/link";
import pathTool from "path";
import { memo } from "react";
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

const Directory = memo(function RenderDirecotry({
  directory,
}: {
  directory: Directory;
}) {
  return (
    <ul>
      {directory.files.map((file) => {
        let href = "/docs/";
        if (directory.slug.length > 0) href += directory.slug.join("/") + "/";
        href += file;
        return (
          <NextLink key={file} href={href}>
            <li>
              <Link>{file}</Link>
            </li>
          </NextLink>
        );
      })}
      {directory.subDirectories.map((directory) => {
        return (
          <div key={directory.dirName}>
            <li>{directory.dirName}</li>
            <Directory directory={directory} />
          </div>
        );
      })}
    </ul>
  );
});

interface Props {
  frontMatter: Record<string, any>;
  mdxSource: MDXRemoteSerializeResult;
  slug: Array<string>;
  filePath: Directory;
}

const PageNavigator = ({ filePath }: { filePath: Directory }) => {
  return <Directory directory={filePath} />;
};

const StyledDocsLayout = styled.div`
  display: flex;
  flex-direction: row;
`;
const DocContent = styled.div`
  margin: 16px;
`;

const components = { Image };
const DocPage = ({ mdxSource, filePath }: Props) => {
  return (
    <StyledDocsLayout>
      <PageNavigator filePath={filePath} />
      <DocContent>
        <MDXRemote {...mdxSource} components={components} />
      </DocContent>
    </StyledDocsLayout>
  );
};

export default DocPage;
