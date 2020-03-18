import React from "react"
import { Flex, Box, Button, Heading, Text } from "theme-ui"

import { api_prefix } from "../globals"

const ExploreData = ({ project, files }) => {
  const matrix_files = files.filter(file => (
    file.file_format === "loom"
  ))

  return (
      <Flex>
        <Button onClick={() => fetch(`${api_prefix}/files/${matrix_files[0].id}`)
            .then(response => response.json()) // parse JSON from request
            .then(resultData => {
              console.log(resultData)
              window.open(resultData.url)
            })}
        >Download matrix</Button>
        <Button>Open in cellxgene</Button>
      </Flex>
  )
}
export default ExploreData
