declare module "census-projects.json" {
  export interface StaticProject
    extends Partial<import("src/common/queries/censusDirectory").Project> {
    notebook_links?: [string, string][];
    tier: number;
  }
  const content: StaticProject[];
  export default content;
}
