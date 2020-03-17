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
          width: "200px", // Narrower so Project Name has more width; this accommodates the longest organ name in Montserrat size 14px
          margin: "0px",
          marginRight: 6,
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
          width: "170px", // Narrower so Project Name has more width; this accommodates the longest assay name in Montserrat size 14px
          margin: "0px",
          marginRight: 6,
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
          width: "120px", // Narrower so Project Name has more width; this accommodates the longest species name in Montserrat size 14px
          margin: "0px",
          marginRight: 6,
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
          width: "100px", // Narrower so Project Name has more width; this accommodates "Unspecified" in Montserrat size 14px
          margin: "0px",
          marginRight: 8,
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
          width: `calc(100% - 830px)`, // calculating the remaining width to allot to project name, after accounting for page margins, and the other columns' widths + their margins
          margin: "0px",
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
