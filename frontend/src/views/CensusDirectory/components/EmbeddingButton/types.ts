import { type StaticProject } from "census-projects.json";
import { type Project } from "src/common/queries/censusDirectory";

export interface EmbeddingButtonProps {
  project: StaticProject | Project;
}
