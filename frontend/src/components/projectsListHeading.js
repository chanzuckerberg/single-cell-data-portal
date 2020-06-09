import React from "react"
import { Flex, Box } from "theme-ui"

const ProjectsListHeading = ({ projects }) => {
  return (
    <Flex sx={{ mb: [4] }}>
      <Box
        sx={{
          fontSize: [2],
          boxSizing: "border-box",
          flexGrow: 0,
          flexShrink: 0,
          width: `16%`, // Narrower so Project Name has more width; account for right margin of "6" in array = 48px
          paddingRight: 6,
          overflow: "hidden", // Or flex might break
          listStyle: "none",
          fontWeight: 700,
        }}
      >
        Organ
      </Box>

      <Box
        sx={{
          fontSize: [2],
          boxSizing: "border-box",
          flexGrow: 0,
          flexShrink: 0,
          width: `16%`, // Narrower so Project Name has more width; account for right margin of "6" in array = 48px
          paddingRight: 6,
          overflow: "hidden", // Or flex might break
          listStyle: "none",
          fontWeight: 700,
        }}
      >
        Assay
      </Box>

      <Box
        sx={{
          fontSize: [2],
          boxSizing: "border-box",
          flexGrow: 0,
          flexShrink: 0,
          width: `16%`, // Narrower so Project Name has more width; account for right margin of "6" in array = 48px
          paddingRight: 6,
          overflow: "hidden", // Or flex might break
          listStyle: "none",
          fontWeight: 700,
        }}
      >
        Species
      </Box>

      <Box
        sx={{
          fontSize: [2],
          boxSizing: "border-box",
          flexGrow: 0,
          flexShrink: 0,
          width: `16%`, // Narrower so Project Name has more width; account for right margin of "6" in array = 48px
          paddingRight: 7,
          overflow: "hidden", // Or flex might break
          listStyle: "none",
          fontWeight: 700,
          textAlign: "right", // since it holds numbers
        }}
      >
        Cells
      </Box>

      <Box
        sx={{
          fontSize: [2],
          boxSizing: "border-box",
          flexGrow: 0,
          flexShrink: 0,
          width: `36%`, // Narrower so Project Name has more width; account for right margin of "6" in array = 48px
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
