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
              paddingBottom: [5],
              paddingTop: "10px",
              borderTop: "#e5e5e5 1px solid",
            }}
            key={project.id}
          >
            <Box
              sx={{
                fontSize: [1],
                boxSizing: "border-box",
                flexGrow: 0,
                flexShrink: 0,
                width: "200px", // Narrower so Project Name has more width
                margin: "0px",
                marginRight: 6,                overflow: "hidden", // Or flex might break
                listStyle: "none",
              }}
            >
              <Flex sx={{ flexDirection: "column" }}>
                {project.tissues.map(tissue=> (
                  <Box key={tissue}>{tissue}</Box>
                ))}
              </Flex>
            </Box>
            <Box
              sx={{
                fontSize: [1],
                boxSizing: "border-box",
                flexGrow: 0,
                flexShrink: 0,
                width: "170px", // Narrower so Project Name has more width
                margin: "0px",
                marginRight: 6,                overflow: "hidden", // Or flex might break
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
                width: "120px", // Narrower so Project Name has more width
                margin: "0px",
                marginRight: 6,                overflow: "hidden", // Or flex might break
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
                width: "100px", // Narrower so Project Name has more width
                margin: "0px",
                marginRight: 8,
                overflow: "hidden", // Or flex might break
                listStyle: "none",
                textAlign: "right", // since it holds numbers
              }}
            >
              {project.cell_count || "Unspecified"}
            </Box>
            <Box
              sx={{
                fontSize: [1],
                boxSizing: "border-box",
                flexGrow: 0,
                flexShrink: 0,
                width: `calc(100% - 830px)`, // calculating the remaining width to allot to project name, after accounting for page margins, and the other columns' widths + their margins
                margin: "0px",
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
