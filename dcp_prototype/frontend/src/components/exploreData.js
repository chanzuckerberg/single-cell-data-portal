import React from "react"
import { Flex, Button, Text } from "theme-ui"

import { api_prefix } from "../globals"

const ExploreData = ({ project, files }) => {

  return (
    !files || !files.length ? <Text>No files available for download</Text> : (
      <Flex>
        <Button onClick={() => fetch(`${api_prefix}/files/${files[0].id}`)
            .then(response => response.json()) // parse JSON from request
            .then(resultData => {
              window.open(resultData.url)
            })}>
            Download matrix
        </Button>
        <Button>Open in cellxgene</Button>
      </Flex>
    )
  )
}
export default ExploreData
