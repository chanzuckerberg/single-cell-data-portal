import fs from "fs";
import matter from "gray-matter";
import { GetStaticPaths } from "next";
import { MDXRemote, MDXRemoteSerializeResult } from "next-mdx-remote";
import { serialize } from "next-mdx-remote/serialize";
import Image from "next/image";
import pathTool from "path";

const DOC_SITE_FOLDER_NAME = "doc-site";

interface Directory {
  dirName: string;
  files: Array<string>;
  subDirectories: Array<Directory>;
}

const CACHED_FILE_PATHS = new Map<Array<string>, Directory>();
function filePaths(...root: Array<string>): Directory {
  const cacheStore = CACHED_FILE_PATHS && CACHED_FILE_PATHS.get(root);
  if (cacheStore) return cacheStore;

  // const paths = fs.readdirSync(pathTool.join(DOC_SITE_FOLDER_NAME, ...root));
  // const builtPaths = paths.reduce((acc, path) => {
  //   const builtFullPath = pathTool.join(DOC_SITE_FOLDER_NAME, ...root, path);
  //   if (fs.lstatSync(builtFullPath).isFile()) {
  //     if (path.endsWith(".mdx")) acc.push([...root, path]);
  //   } else {
  //     acc.push(...filePaths(...root, path));
  //   }
  //   return acc;
  // }, {} as Directory);};
  const newDirectory = {
    dirName: root[root.length - 1] ?? "",
    files: [],
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
  CACHED_FILE_PATHS.set(root, newDirectory);
  return newDirectory;
}

// need to fix the index path
// There are incorrect paths that are reminants from previous recursive calls
function generatePaths(
  directory: Directory,
  slugs = new Array<{ params: { slug: Array<string> } }>(),
  indexPath = new Array<string>()
): Array<{ params: { slug: Array<string> } }> {
  const newIndexPath = [...indexPath];
  if (directory.dirName) newIndexPath.push(directory.dirName);
  directory.files.forEach((file) => {
    slugs.push({ params: { slug: [...newIndexPath, file] } });
  });
  if (directory.subDirectories.length > 0) {
    directory.subDirectories.forEach((subDirectory) => {
      slugs.push(...generatePaths(subDirectory, slugs, newIndexPath));
    });
  }
  return slugs;
}

export const getStaticPaths: GetStaticPaths = async () => {
  const filePath = filePaths();
  // console.log(filePath);
  const paths = generatePaths(filePath);
  paths.forEach((path) => {
    console.log(path.params.slug);
  });

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
  const { data: frontMatter, content } = matter(markdownWithMeta);
  const mdxSource = await serialize(content);
  return {
    props: {
      frontMatter,
      mdxSource,
      slug,
    },
  };
};

const components = { Image };

interface Props {
  frontMatter: Record<string, any>;
  mdxSource: MDXRemoteSerializeResult;
  slug: Array<string>;
}

const PageNavigator = ({ pages }) => {
  return <div></div>;
};

const BlogPage = ({ mdxSource }: Props) => {
  return (
    <>
      <MDXRemote {...mdxSource} components={components} />
    </>
  );
};

export default BlogPage;
