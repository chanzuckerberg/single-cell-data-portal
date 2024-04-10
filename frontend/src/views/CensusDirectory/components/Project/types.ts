import { type StaticProject } from "census-projects.json";
import { type Project } from "src/common/queries/censusDirectory";
import { type ClobberedProjects } from "../../utils";

export type UnionProject = StaticProject | Project;
export interface ProjectProps {
  clobberedProjects: ClobberedProjects[number];
}
