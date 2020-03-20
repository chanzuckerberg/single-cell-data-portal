import React from "react"
import { Flex, Box, Heading, Text } from "theme-ui"

const ProjectOverview = ({ project }) => {
  return (
    <Box>
      <Heading as="h3" sx={{ mb: 4 }}>
        {project.title}
      </Heading>
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
          <Box>
            <Heading as="h6" sx={{ mb: [2] }}>
              Publications
            </Heading>
            <Text>{project.publication_title}</Text>
          </Box>
          <Box>
            <Heading as="h6" sx={{ mb: [2] }}>
              Contributors
            </Heading>
            <Text>{project.contributors[0].first_name}</Text>
          </Box>
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
