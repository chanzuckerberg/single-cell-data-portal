import React from "react"
import { Flex, Button, Text } from "theme-ui"

import { api_prefix } from "../globals"

const ExploreData = ({ project, files }) => {
  const cxg_url = "http://cellxgene-dev.us-west-2.elasticbeanstalk.com/"

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
        <Button onClick={() => window.open(cxg_url + project.title + ".cxg")}>
            Visualize in cellxgene
        </Button>
      </Flex>
    )
  )
}
export default ExploreData
