import { type StaticProject } from "census-projects.json";
import {
  type Project,
  type ProjectResponse,
} from "src/common/queries/censusDirectory";

export interface ProjectProps {
  project: StaticProject | Project;
  id?: keyof ProjectResponse;
}
