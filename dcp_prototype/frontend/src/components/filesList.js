import React from "react"
import { Heading, Text, Flex, Box } from "theme-ui"

const FilesList = ({ files }) => {
  console.log(files)
  return (
    <Box>
      {files.map(file => {
        return (
          <Flex key={file.id} sx={{ justifyContent: "space-between" }}>
            <Text sx={{ fontSize: 1 }}>Filename: {file.filename}</Text>
            <Text sx={{ fontSize: 1 }}>{file.species}</Text>
            <Text sx={{ fontSize: 1 }}>{file.file_size}</Text>
            <Text sx={{ fontSize: 1 }}>
              {file.library_construction_method_ontology}
            </Text>
            <Text sx={{ fontSize: 1 }}>{file.tissue_ontology}</Text>
            <Text sx={{ fontSize: 1 }}>{file.file_format}</Text>
          </Flex>
        )
      })}
    </Box>
  )
}

export default FilesList

// id: "HCA-AnalysisFile-001a3fdd-c9e6-4871-9e48-a840d52ecdf9"
// filename: "b29eeb85-53ff-4786-9101-3e241e6dc250_rsem.bam"
// file_format: "bam"
// file_size: 0
// species: "homo sapiens"
// library_construction_method_ontology: ""
// tissue_ontology: ""
