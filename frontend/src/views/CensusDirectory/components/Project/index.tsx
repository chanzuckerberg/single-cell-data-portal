import Link from "next/link";
import { type Project } from "src/common/queries/censusDirectory";
import {
  ProjectContainer,
  ProjectDetails,
  ProjectTitle,
  ProjectSubmitter,
  ProjectDesctiption,
  DetailsContainer,
  ProjectButtons,
  StyledButton,
} from "../../styles";
import DetailItem from "../DetailItem";

import EmbeddingButton from "../EmbeddingButton";

import { ProjectProps } from "./types";
import { useConnect } from "./connect";

const Project = ({ project, id }: ProjectProps) => {
  const { date, projectNotebookLinks } = useConnect({ project, id });

  return (
    <ProjectContainer key={project.title}>
      <ProjectDetails>
        <ProjectTitle>{project.title}</ProjectTitle>
        <ProjectSubmitter>{project.contact_affiliation}</ProjectSubmitter>
        <ProjectDesctiption>{project.description}</ProjectDesctiption>
        <DetailsContainer>
          <DetailItem label="contact" link={project.contact_email}>
            {project.contact_name}
          </DetailItem>
          <DetailItem label="publication" link={project.publication_link}>
            {project.publication_info}
          </DetailItem>
          <DetailItem label="Last Updated">{date}</DetailItem>
        </DetailsContainer>
        <DetailsContainer>
          <DetailItem label="Census Version">
            {project.census_version}
          </DetailItem>
          <DetailItem label="experiment">{project.experiment_name}</DetailItem>
          <DetailItem label="measurement">
            {project.measurement_name}
          </DetailItem>
          <DetailItem label="embedding">{project.data_type}</DetailItem>
          {projectNotebookLinks?.map((link) => (
            <DetailItem label="notebook" link={link[1]} key={link[1]}>
              {link[0]}
            </DetailItem>
          ))}
        </DetailsContainer>
      </ProjectDetails>
      <ProjectButtons>
        <EmbeddingButton />
        {project.model_link && (
          <Link href={project.model_link}>
            <StyledButton sdsType="primary" sdsStyle="square">
              Model
            </StyledButton>
          </Link>
        )}
      </ProjectButtons>
    </ProjectContainer>
  );
};
export default Project;
