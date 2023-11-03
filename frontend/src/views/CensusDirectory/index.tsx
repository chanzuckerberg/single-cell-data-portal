import {
  Button,
  fontBodyS,
  fontCapsXxxs,
  fontHeaderL,
  fontHeaderXl,
  fontHeaderXxl,
} from "@czi-sds/components";
import styled from "@emotion/styled";
import Link from "next/link";
import React from "react";
import {
  Project,
  ProjectResponse,
  useProjects,
} from "src/common/queries/censusDirectory";
import {
  fontWeightBold,
  fontWeightRegular,
  fontWeightSemibold,
  gray400,
  spacesDefault,
  spacesL,
  spacesXl,
  spacesXxs,
  textSecondary,
} from "src/common/theme";

import staticProjects, { type StaticProject } from "census-projects.json";
import notebookLinks from "census-notebook-links.json";
import { useRouter } from "next/router";

function DetailItem(props: { label: string; children: string; link?: string }) {
  const ItemContainer = styled.div`
    display: flex;
    flex-direction: column;
    gap: ${spacesXxs}px;
  `;

  const ItemLabel = styled.div`
    ${fontCapsXxxs}
    font-weight: ${fontWeightSemibold};
    font-feature-settings:
      "clig" off,
      "liga" off;
    color: ${gray400};
  `;

  return (
    <ItemContainer>
      <ItemLabel>{props.label}</ItemLabel>
      {props.link ? (
        <Link href={props.link} passHref>
          {props.children}
        </Link>
      ) : (
        props.children
      )}
    </ItemContainer>
  );
}

function CensusDirectory() {
  const Content = styled.div`
    display: flex;
    flex-direction: column;
    margin: 80px auto;
    max-width: 1200px;
  `;

  const Header = styled.h1`
    ${fontHeaderXxl}
    margin-bottom: ${spacesDefault}px;
    font-weight: ${fontWeightBold};
  `;

  const Paragraph = styled.p`
    ${fontBodyS}
    font-weight: ${fontWeightRegular};
    margin-bottom: 0;
  `;

  const DirectoryDescription = styled(Paragraph)`
    margin-bottom: 80px;
  `;

  const TierContainer = styled.div`
    margin-bottom: 120px;
  `;

  const TierTitle = styled.h3`
    ${fontHeaderXl}
    margin-bottom: ${spacesDefault}px;
    font-weight: ${fontWeightSemibold};
  `;

  const TierDescription = styled.p`
    ${fontBodyS}
    color: ${textSecondary};
    font-weight: ${fontWeightRegular};
    margin-bottom: 0;
  `;

  const ProjectTitle = styled.h4`
    ${fontHeaderL}
    font-weight: ${fontWeightSemibold};
    margin-bottom: ${spacesDefault}px;
  `;

  const ProjectSubmitter = styled.h4`
    ${fontBodyS}
    font-weight: ${fontWeightSemibold};
    margin-bottom: ${spacesDefault}px;
  `;

  const ProjectDesctiption = styled(Paragraph)`
    max-width: 85ch;
  `;

  const ProjectContainer = styled.div`
    display: flex;
    flex-direction: row;
    justify-content: space-between;
    margin-top: 20px;
  `;
  const ProjectButtons = styled.div`
    display: flex;
    flex-direction: row;
    gap: ${spacesDefault}px;
  `;
  const ProjectDetails = styled.div`
    display: flex;
    flex-direction: column;
  `;
  const DetailsContainer = styled.div`
    display: flex;
    flex-direction: row;
    gap: ${spacesXl}px;
    margin-top: ${spacesL}px;
  `;

  const StyledButton = styled(Button)`
    font-weight: ${fontWeightSemibold};
    min-width: 80px;
  `;

  const { data: projects } = useProjects();

  const hostedProjects = Object.entries(projects || {}).filter(
    ([_, project]) => !project.revised_by
  );

  const communityProjects = Object.values(staticProjects).filter(
    (project) => project.tier === 1
  );
  const maintainedProjects = Object.values(staticProjects).filter(
    (project) => project.tier === 3
  );

  const router = useRouter();

  const renderDetailItem = (label: string, value?: string, link?: string) => {
    return value ? (
      <DetailItem label={label} link={link}>
        {value}
      </DetailItem>
    ) : null;
  };

  const renderProject = (
    project: StaticProject | Project,
    id?: keyof ProjectResponse
  ) => (
    <ProjectContainer key={project.title}>
      <ProjectDetails>
        <ProjectTitle>{project.title}</ProjectTitle>
        <ProjectSubmitter>{project.contact_affiliation}</ProjectSubmitter>
        <ProjectDesctiption>{project.description}</ProjectDesctiption>
        <DetailsContainer>
          {renderDetailItem(
            "contact",
            project.contact_name,
            project.contact_email
          )}
          {renderDetailItem(
            "publication",
            project.publication_info,
            project.publication_link
          )}
          {renderDetailItem(
            "Last Updated",
            //convert date to month, day, year
            new Date(
              project.last_updated || project.submission_date || ""
            ).toLocaleDateString("en-US", {
              dateStyle: "long",
            })
          )}
        </DetailsContainer>
        <DetailsContainer>
          {renderDetailItem("Census Version", project.census_version)}
          {renderDetailItem("experiment", project.experiment_name)}
          {renderDetailItem("measurement", project.measurement_name)}
          {renderDetailItem("embedding", project.data_type)}
          {"notebook_links" in project
            ? project.notebook_links?.map((link) =>
                renderDetailItem("notebook", link[0], link[1])
              )
            : id &&
              notebookLinks[id]?.map((link) =>
                renderDetailItem("notebook", link[0], link[1])
              )}
        </DetailsContainer>
      </ProjectDetails>
      <ProjectButtons>
        <StyledButton sdsType="secondary" sdsStyle="square">
          Embedding
        </StyledButton>
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
  return (
    <Content>
      <Header>Census Directory</Header>
      <DirectoryDescription>
        Habitant sit tristique pharetra at. Quis ipsum morbi pharetra venenatis
        amet purus aliquam nunc. Mi feugiat elementum nec sagittis enim. Turpis
        purus mi tellus leo in vestibulum enim varius. Ut dictum lobortis in
        non. Sed rhoncus enim pharetra pulvinar semper faucibus ut at sapien.
        Parturient pharetra amet sit facilisis sagittis quis. Dignissim
        fermentum consectetur fames vulputate semper neque est non pharetra.
        Amet et elementum neque turpis hac bibendum ac id ipsum.
      </DirectoryDescription>
      {maintainedProjects.length > 0 && (
        <TierContainer>
          <TierTitle>Census Partners</TierTitle>
          <TierDescription>
            Ut nisi non lorem adipiscing. Orci tellus quisque quam ac purus
            vitae. Aliquet quis egestas viverra nulla quis lectus adipiscing.
          </TierDescription>
          {maintainedProjects.map((project) => renderProject(project))}
        </TierContainer>
      )}
      {hostedProjects.length > 0 && (
        <TierContainer>
          <TierTitle>Tier 2</TierTitle>
          <TierDescription>
            Ut nisi non lorem adipiscing. Orci tellus quisque quam ac purus
            vitae. Aliquet quis egestas viverra nulla quis lectus adipiscing.
          </TierDescription>
          {hostedProjects.map(([id, project]) => renderProject(project, id))}
        </TierContainer>
      )}
      {communityProjects.length > 0 && (
        <TierContainer>
          <TierTitle>Tier 3</TierTitle>
          <TierDescription>
            Ut nisi non lorem adipiscing. Orci tellus quisque quam ac purus
            vitae. Aliquet quis egestas viverra nulla quis lectus adipiscing.
          </TierDescription>
          {communityProjects.map((project) => renderProject(project))}
        </TierContainer>
      )}
    </Content>
  );
}

export default CensusDirectory;
