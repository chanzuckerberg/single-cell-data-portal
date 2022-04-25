import SidebarLayout from "src/components/Layout/components/sidebarLayout";

// export const getStaticProps = async () => {
//   const files = fs.readdirSync(path.join("doc-site"));
//   const docPages = files.reduce((acc, filename) => {
//     if (fs.lstatSync(path.join("doc-site", filename)).isDirectory()) {
//       console.log(path.join("doc-site", filename), " IS DIRECTORY");
//       return acc;
//     } else {
//       console.log(path.join("doc-site", filename), " IS FILE");
//     }
//     const markdownWithMeta = fs.readFileSync(
//       path.join("doc-site", filename),
//       "utf-8"
//     );
//     const { data: frontMatter } = matter(markdownWithMeta);
//     acc.push({
//       frontMatter,
//       slug: filename.split(".")[0],
//     });
//     return acc;
//   }, new Array<{ frontMatter: Record<string, any>; slug: string }>());

Page.Layout = SidebarLayout;

export default Page;
