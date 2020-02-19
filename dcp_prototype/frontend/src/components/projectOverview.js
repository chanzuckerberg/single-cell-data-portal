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
              Contact
            </Heading>
            <Text>{project.contact_name}</Text>
            <Text>{project.contact_institution}</Text>
            <Text>{project.contact_email}</Text>
          </Box>
        </Box>
        <Flex>
          <Box>
            <Heading as="h6" sx={{ mb: [2] }}>
              Project Details
            </Heading>
            <Text>Project Label: {project.label}</Text>
            <Text>Species: {project.species.map(species => species)}</Text>
            <Text>Sample type: {project.sample_type}</Text>
            <Text>Organ Part: {project.organ_part}</Text>
            <Text>
              Analysis Protocols:
              {project.analysis_protocol.map(protocol => protocol)}
            </Text>
            <Text>Cell count: {project.cell_count}</Text>
            <Text>Donor count: {project.donor_count}</Text>
          </Box>
        </Flex>
      </Flex>
    </Box>
  )
}
export default ProjectOverview
