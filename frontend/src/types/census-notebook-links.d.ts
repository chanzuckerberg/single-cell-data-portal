declare module "census-notebook-links.json" {
  const content: {
    [accessionID: string]: [title: string, link: string][];
  };
  export default content;
}
