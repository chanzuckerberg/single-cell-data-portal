import React, { FC } from "react";
import { Link } from "gatsby";
import { Flex, Box } from "theme-ui";
import ProjectsListHeading from "./projectsListHeading";
import { Project } from "../common/entities";
import { SystemStyleObject } from "@styled-system/css";

interface Props {
  projects: Project[];
}

const boxBaseStyle: SystemStyleObject = {
  fontSize: [1],
  boxSizing: "border-box",
  flexGrow: 0,
  flexShrink: 0,
  width: `16%`, // Narrower so Project Name has more width; account for right margin of "6" in array = 48px
  paddingRight: 6,
  overflow: "hidden", // Or flex might break
  listStyle: "none",
};

const StyledBox: FC = props => {
  return <Box {...props} sx={boxBaseStyle} />;
};

const CellCountBox: FC = props => {
  const style: SystemStyleObject = {
    ...boxBaseStyle,
    paddingRight: 7,
    textAlign: "right", // since it holds numbers
    paddingTop: "10px",
    lineHeight: "1.5",
  };

  return <Box {...props} sx={style} />;
};

const ProjectLinkBox: FC = props => {
  const style: SystemStyleObject = {
    ...boxBaseStyle,
    width: `36%`, // Narrower so Project Name has more width; account for right margin of "6" in array = 48px
    paddingRight: "unset",
    paddingTop: "10px",
    lineHeight: "1.5",
  };

  return <Box {...props} sx={style} />;
};

const ProjectDetailBox: FC = props => {
  return (
    <Box
      {...props}
      sx={{
        paddingTop: "10px",
        lineHeight: "1.5",
      }}
    />
  );
};

interface ProjectDetailProps {
  entities: string[];
}

const ProjectDetail: FC<ProjectDetailProps> = ({ entities }) => {
  return (
    <StyledBox>
      <Flex sx={{ flexDirection: "column" }}>
        {entities.map(entity => (
          <ProjectDetailBox key={entity}>{entity}</ProjectDetailBox>
        ))}
      </Flex>
    </StyledBox>
  );
};

interface ProjectProps {
  project: Project;
}

const ProjectComponent: FC<ProjectProps> = ({ project }) => {
  return (
    <Flex
      sx={{
        fontSize: [1],
        paddingBottom: [5],
        borderTop: "#e5e5e5 1px solid",
      }}
      key={project.id}
    >
      <StyledBox>
        <Flex sx={{ flexDirection: "column" }}>
          {project.organs.map(organ => (
            <ProjectDetailBox key={organ}>
              {organ.charAt(0).toUpperCase() + organ.substring(1)}
            </ProjectDetailBox>
          ))}
        </Flex>
      </StyledBox>
      <ProjectDetail entities={project.assays} />
      <ProjectDetail entities={project.species} />
      <CellCountBox>{project.cell_count || "Unspecified"}</CellCountBox>
      <ProjectLinkBox>
        <Link to={`/project/?id=${project.id}`}>{project.title}</Link>
      </ProjectLinkBox>
    </Flex>
  );
};

const ProjectsList: FC<Props> = ({ projects }) => {
  return (
    <>
      <ProjectsListHeading />
      {projects.map(project => (
        <ProjectComponent key={project.id} project={project} />
      ))}
    </>
  );
};

export default ProjectsList;
