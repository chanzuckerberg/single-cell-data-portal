import React from "react"
import { Flex, Box } from "theme-ui"

const ProjectsListHeading = ({ projects }) => {
  return (
    <Flex sx={{ mb: [4] }}>
      <Box
        sx={{
          fontSize: [3],
          boxSizing: "border-box",
          flexGrow: 0,
          flexShrink: 0,
          width: "20%", // Default to full width
          // padding: [0],
          overflow: "hidden", // Or flex might break
          listStyle: "none",
          fontWeight: 700,
        }}
      >
        Organ
      </Box>

      <Box
        sx={{
          fontSize: [3],
          boxSizing: "border-box",
          flexGrow: 0,
          flexShrink: 0,
          width: "20%", // Default to full width
          // padding: [0],
          overflow: "hidden", // Or flex might break
          listStyle: "none",
          fontWeight: 700,
        }}
      >
        Assay
      </Box>

      <Box
        sx={{
          fontSize: [3],
          boxSizing: "border-box",
          flexGrow: 0,
          flexShrink: 0,
          width: "20%", // Default to full width
          // padding: [0],
          overflow: "hidden", // Or flex might break
          listStyle: "none",
          fontWeight: 700,
        }}
      >
        Species
      </Box>

      <Box
        sx={{
          fontSize: [3],
          boxSizing: "border-box",
          flexGrow: 0,
          flexShrink: 0,
          width: "20%", // Default to full width
          // padding: [0],
          overflow: "hidden", // Or flex might break
          listStyle: "none",
          fontWeight: 700,
        }}
      >
        Cells
      </Box>

      <Box
        sx={{
          fontSize: [3],
          boxSizing: "border-box",
          flexGrow: 0,
          flexShrink: 0,
          width: "40%", // Default to full width
          // padding: [0],
          overflow: "hidden", // Or flex might break
          listStyle: "none",
          fontWeight: 700,
        }}
      >
        Project name
      </Box>
    </Flex>
  )
}

export default ProjectsListHeading
