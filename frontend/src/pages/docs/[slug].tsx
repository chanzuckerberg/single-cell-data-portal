import { H1 } from "@blueprintjs/core";
import fs from "fs";
import matter from "gray-matter";
import { GetStaticPaths } from "next";
import { MDXRemote } from "next-mdx-remote";
import { MDXRemoteSerializeResult } from "next-mdx-remote/dist/types";
import { serialize } from "next-mdx-remote/serialize";
import path from "path";

export const getStaticPaths: GetStaticPaths = async () => {
  const files = fs.readdirSync(path.join("doc-site"));
  const paths = files.map((filename) => ({
    params: {
      slug: filename.replace(".mdx", ""),
    },
  }));
  return {
    fallback: false,
    paths,
  };
};

export const getStaticProps = async ({
  params: { slug },
}: {
  params: { slug: string };
}): Promise<{
  props: Props;
}> => {
  const markdownWithMeta = fs.readFileSync(
    path.join("doc-site", slug + ".mdx"),
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

interface Props {
  frontMatter: Record<string, any>;
  mdxSource: MDXRemoteSerializeResult;
  slug: string;
}

const BlogPage = ({ frontMatter, mdxSource }: Props) => {
  return (
    <>
      <H1>THIS HEADER WRAPS EVERY SINGLE DOC PAGE</H1>
      <sub>{frontMatter.toString()}</sub>
      <MDXRemote {...mdxSource} />
    </>
  );
};

export default BlogPage;
