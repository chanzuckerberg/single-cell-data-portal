import { ClobberedProjects } from "src/views/CensusDirectory/utils";
import { UnionProject } from "../../../types";

export interface EmbeddingButtonProps {
  project: UnionProject;
  uniqueMetadata: ClobberedProjects[number][0];
}
