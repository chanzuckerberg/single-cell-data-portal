import React from "react"
import { Flex, Box, Heading, Text } from "theme-ui"
import ExploreData from "./exploreData"

const ProjectOverview = ({ project, files, isAuthenticated }) => {
  const publications = !project.publication_title ? null : (
      <Box sx={{mb: [4]}}>
        <Heading as="h6" sx={{ mb: [2] }}>
          Publications
        </Heading>
        <Text>{project.publication_title}</Text>
      </Box>
  )
  const contributors = !project.contributors.length ? null : (
      <Flex>
        <Box>
          <Heading as="h6" sx={{ mb: [2] }}>
            Contributors
          </Heading>
          {project.contributors.map(c => <Text key={c.name}>{c.name}</Text>)}
        </Box>
        <Box sx={{ml: [2]}}>
          <Heading as="h6" sx={{ mb: [2] }}>
            Institutions
          </Heading>
          <Text>{project.contributors[0].institution}</Text>
        </Box>
      </Flex>
  )
  return (
    <Box>
      <Heading as="h3" sx={{ mb: 4 }}>
        {project.title}
      </Heading>
      {!project ? null : (
        <Flex sx={{flexDirection: "column"}}>
          <ExploreData project={project} files={files} isAuthenticated={isAuthenticated}/>
        </Flex>
      )}
      <Heading as="h3" sx={{ mb: 4 }}>
        Project Information
      </Heading>
      <Flex
        sx={{
          fontSize: [1],
          mb: [4],
        }}
      >
        <Box sx={{ maxWidth: "30em", mr: 5 }}>
          <Box sx={{ mb: [4] }}>
            <Heading as="h6" sx={{ mb: [2] }}>
              Description
            </Heading>
            <Text>{project.description}</Text>
          </Box>
          {publications}
          {contributors}
        </Box>
        <Flex>
          <Box>
            <Heading as="h6" sx={{ mb: [2] }}>
              Project Details
            </Heading>
            <Text>Species: {project.species.join(", ")}</Text>
            <Text>Organs: {project.organs.join(", ")}</Text>
            <Text>Sample Categories: {project.biosample_categories.join(", ")}</Text>
            <Text>Disease Status: {project.diseases.join(", ")}</Text>
            <Text>Library Construction methods: {project.assays.join(", ")}</Text>
            <Text>Paired end: {project.paired_end.join(", ")}</Text>
            <Text>Cell count estimate: {project.cell_count}</Text>
          </Box>
        </Flex>
      </Flex>
    </Box>
  )
}
export default ProjectOverview
