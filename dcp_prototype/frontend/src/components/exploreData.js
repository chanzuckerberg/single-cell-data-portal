import React from "react"
import { Flex, Box, Text } from "theme-ui"
import StyledButton from "./styledButton"
import LoginSignup from "./login-signup"

import { api_prefix } from "../globals"

const ExploreData = ({ project, files, isAuthenticated }) => {
  const cxg_url = "http://cellxgene-dev.us-west-2.elasticbeanstalk.com/"
  const disabled = !((files && files.length && isAuthenticated) === true)

  return (
  <Box sx={{mb: 2}}>
    <Flex sx={{mb: 1}}>
      <StyledButton label="Download matrix"
                    disabled={disabled}
                    onclick={() => fetch(`${api_prefix}/files/${files[0].id}`)
                                    .then(response => response.json())
                                    .then(resultData => {
                                      window.open(resultData.url)
                                    })}>
        Download matrix
      </StyledButton>
      <StyledButton label="Visualize in cellxgene"
                    disabled={disabled}
                    onclick={() => window.open(cxg_url + project.title + ".cxg")}>
          Visualize in cellxgene
      </StyledButton>
    </Flex>
    { isAuthenticated ? null : (
      <Box>
        <LoginSignup/>
      </Box>
    )}
    { files && !files.length && isAuthenticated ? (
      <Text sx={{color: "primary", fontSize: 1}}>
        No files available for download
      </Text>
    ) : null }
  </Box>
  )
}
export default ExploreData
