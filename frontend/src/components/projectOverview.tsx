import React, { FC } from "react";
import { Box, Flex, Heading, Text } from "theme-ui";
import { File, Project } from "../common/entities";
import ExploreData from "./exploreData";

interface Props {
  project: Project;
  files: [File];
  isAuthenticated: boolean;
}

const StyledHeading: FC = props => {
  return <Heading {...props} as="h6" sx={{ mb: [2] }} />;
};

interface ProjectDetailProps {
  heading: string;
  text: string;
}

const ProjectDetail: FC<ProjectDetailProps> = ({ heading, text }) => {
  return (
    <Box sx={{ mb: [4] }}>
      <StyledHeading>{heading}</StyledHeading>
      <Text>{text}</Text>
    </Box>
  );
};

const ProjectOverview: FC<Props> = ({ project, files, isAuthenticated }) => {
  const description = (
    <ProjectDetail heading="Description" text={project.description} />
  );

  const publications = !project.publication_title ? null : (
    <ProjectDetail heading="Publications" text={project.publication_title} />
  );

  const contributors = !project.contributors.length ? null : (
    <Flex>
      <Box
        sx={{
          width: "50%",
        }}
      >
        <StyledHeading>Contributors</StyledHeading>
        {project.contributors.map(c => (
          <Text key={c.name}>{c.name}</Text>
        ))}
      </Box>
      <Box
        sx={{
          ml: [2],
          marginLeft: [6],
        }}
      >
        <StyledHeading>Institutions</StyledHeading>
        <Text>{project.contributors[0].institution}</Text>
      </Box>
    </Flex>
  );

  return (
    <Box>
      <Heading
        as="h3"
        sx={{
          mb: 4,
          maxWidth: 0,
        }}
      >
        {project.title}
      </Heading>
      {!project ? null : (
        <Flex sx={{ flexDirection: "column" }}>
          <ExploreData
            project={project}
            files={files}
            isAuthenticated={isAuthenticated}
          />
        </Flex>
      )}
      <Heading as="h3" sx={{ mb: 4 }}>
        Project Information
      </Heading>
      <Flex
        sx={{
          fontSize: [1],
          mb: [4],
          maxWidth: 0,
        }}
      >
        <Box
          sx={{
            maxWidth: "60%",
            mr: 5,
          }}
        >
          {description}
          {publications}
          {contributors}
        </Box>
        <Flex>
          <Box>
            <StyledHeading>Project Details</StyledHeading>
            <Text>
              <strong>Species: </strong>
              {project.species.join(", ")}
            </Text>
            <Text>
              <strong>Organs: </strong>
              {project.organs.join(", ")}
            </Text>
            <Text>
              <strong>Sample Categories: </strong>
              {project.biosample_categories.join(", ")}
            </Text>
            <Text>
              <strong>Disease Status: </strong>
              {project.diseases.join(", ")}
            </Text>
            <Text>
              <strong>Library Construction methods: </strong>
              {project.assays.join(", ")}
            </Text>
            <Text>
              <strong>Paired end: </strong>
              {project.paired_end.join(", ")}
            </Text>
            <Text>
              <strong>Cell count estimate: </strong>
              {project.cell_count}
            </Text>
          </Box>
        </Flex>
      </Flex>
    </Box>
  );
};
export default ProjectOverview;
