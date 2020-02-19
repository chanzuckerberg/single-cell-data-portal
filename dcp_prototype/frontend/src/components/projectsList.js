import { Link } from "gatsby"
import React from "react"
import { Flex, Box } from "theme-ui"
import ProjectsListHeading from "./projectsListHeading"

const ProjectsList = ({ projects }) => {
  return (
    <>
      <ProjectsListHeading />
      {projects.map(project => {
        return (
          <Flex
            sx={{
              fontSize: [1],
              mb: [4],
            }}
            key={project.id}
          >
            <Box
              sx={{
                fontSize: [1],
                boxSizing: "border-box",
                flexGrow: 0,
                flexShrink: 0,
                width: "20%", // Default to full width
                // padding: [0],
                overflow: "hidden", // Or flex might break
                listStyle: "none",
              }}
            >
              <Flex sx={{ flexDirection: "column" }}>
                {project.assays.map(assay => (
                  <Box key={assay}>{assay}</Box>
                ))}
              </Flex>
            </Box>
            <Box
              sx={{
                fontSize: [1],
                boxSizing: "border-box",
                flexGrow: 0,
                flexShrink: 0,
                width: "20%", // Default to full width
                // padding: [0],
                overflow: "hidden", // Or flex might break
                listStyle: "none",
              }}
            >
              <Flex sx={{ flexDirection: "column" }}>
                {project.species.map(species => (
                  <Box key={species}>{species}</Box>
                ))}
              </Flex>
            </Box>
            <Box
              sx={{
                fontSize: [1],
                boxSizing: "border-box",
                flexGrow: 0,
                flexShrink: 0,
                width: "20%", // Default to full width
                // padding: [0],
                overflow: "hidden", // Or flex might break
                listStyle: "none",
              }}
            >
              {project.total_cells || "unknown"}
            </Box>
            <Box
              sx={{
                fontSize: [1],
                boxSizing: "border-box",
                flexGrow: 0,
                flexShrink: 0,
                width: "40%", // Default to full width
                // padding: [0],
                overflow: "hidden", // Or flex might break
                listStyle: "none",
              }}
            >
              <Link to={`/project?id=${project.id}`}>{project.title}</Link>
            </Box>
          </Flex>
        )
      })}
    </>
  )
}

export default ProjectsList
