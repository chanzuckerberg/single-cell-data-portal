import fs from "fs";
import matter from "gray-matter";
import { GetStaticPaths } from "next";
import { MDXRemote } from "next-mdx-remote";
import { MDXRemoteSerializeResult } from "next-mdx-remote/dist/types";
import { serialize } from "next-mdx-remote/serialize";
import Image from "next/image";
import pathTool from "path";

const DOC_SITE_FOLDER_NAME = "doc-site";

function recursiveFileGrab(...root: Array<string>): Array<Array<string>> {
  const paths = fs.readdirSync(pathTool.join(DOC_SITE_FOLDER_NAME, ...root));
  return paths.reduce((acc, path) => {
    const builtFullPath = pathTool.join(DOC_SITE_FOLDER_NAME, ...root, path);
    if (fs.lstatSync(builtFullPath).isFile()) {
      acc.push([...root, path]);
    } else {
      acc.push(...recursiveFileGrab(...root, path));
    }
    return acc;
  }, new Array<Array<string>>());
}

export const getStaticPaths: GetStaticPaths = async () => {
  const filePath = recursiveFileGrab();
  const paths = filePath.map((filePath) => {
    const slug = filePath;
    slug[slug.length - 1] = slug[slug.length - 1].replace(".mdx", "");
    return {
      params: {
        slug,
      },
    };
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
  frontMatter: Record<string, unknown>;
  mdxSource: MDXRemoteSerializeResult;
  slug: Array<string>;
}

const BlogPage = ({ mdxSource }: Props): JSX.Element => {
  return (
    <>
      <MDXRemote {...mdxSource} components={components} />
    </>
  );
};

export default BlogPage;
