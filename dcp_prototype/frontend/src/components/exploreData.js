import React from "react"
import { Flex, Box, Text } from "theme-ui"
import StyledButton from "./styledButton"
import LoginSignup from "./login-signup"

import { api_url, cxg_url } from "../globals"
import { useAuth0 } from "../contexts/auth0Context";

const ExploreData = ({ project, files, isAuthenticated }) => {
  const matrixAvailable = (files && files.length && isAuthenticated) === true
  const cxgAvailable = project.cxg_enabled
    const { getTokenSilently } = useAuth0()

  return (
  <Box sx={{mb: 2}}>
    <Flex sx={{mb: 1}}>
      <StyledButton label="Download matrix"
                    disabled={!matrixAvailable}
                    onclick={() => getTokenSilently().then(accessToken =>
                        fetch(`${api_url}/files/${files[0].id}`,{
                            headers: {
                                Authorization: 'Bearer ' + accessToken
                            }
                        }))
                        .then(response => response.json())
                        .then(resultData => {
                            window.open(resultData.url)
                        })}>
        Download matrix
      </StyledButton>
      <StyledButton label="Visualize in cellxgene"
                    disabled={!cxgAvailable}
                    onclick={() => window.open(`${cxg_url}/${project.label}.cxg`)}>
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
