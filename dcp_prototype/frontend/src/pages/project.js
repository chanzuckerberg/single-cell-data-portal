import React, { useState, useEffect } from "react"
import { useAuth0 } from "../contexts/auth0Context"
import Layout from "../components/layout"
import SEO from "../components/seo"
import searchStringAsObj from "../util/searchStringAsObj"
import ProjectOverview from "../components/projectOverview"
import { api_prefix } from "../globals"
import { Flex, Box, Heading } from "theme-ui"

const SecondPage = props => {
  const [project, setProject] = useState(null)
  const [files, setFiles] = useState(null)
  const searchObj = searchStringAsObj(props.location.search.slice(1))
  const id = searchObj?.id
  const { loading, loginWithRedirect, isAuthenticated} = useAuth0()

  useEffect(() => {
    if (id) {
      // TODO send authorization token in requests to backend
      fetch(`${api_prefix}/projects/${id}`)
        .then(response => response.json()) // parse JSON from request
        .then(resultData => {
          setProject(resultData)
        })
    } else {
      return null
    }
  }, [id]) /* don't rerender if this hasn't changed */

  useEffect(() => {
    if (id) {
      // TODO send authorization token in requests to backend
      fetch(`${api_prefix}/projects/${id}/files`)
        .then(response => response.json()) // parse JSON from request
        .then(resultData => {
          setFiles(resultData)
        })
    } else {
      return null
    }
  }, [id]) /* don't rerender if this hasn't changed */

  if(loading){
    return <p>Loading...</p>
  }

  if (!isAuthenticated) {
    loginWithRedirect()
    return <p>Redirecting to login...</p>
  }

  /*
    /projects/id
    /projects/id/files
    /files/id
  */

  return (
    <Layout>
      <SEO title="Projects"/>
      <Heading as="h1" sx={{ mb: 4 }}>
        Explore Project
      </Heading>
      <Flex>
        {!project ? (
          "Loading project..."
        ) : (
          <ProjectOverview project={project}/>
        )}
      </Flex>
      {!project ? null : (
        <Box>
          <Heading as="h3" sx={{ mb: 4 }}>
            Download files
          </Heading>
          <Heading as="h3" sx={{ mb: 4 }}>
            ...
          </Heading>
        </Box>
      )}
    </Layout>
  )
}

export default SecondPage

// <Flex>{!files ? null : <FilesList files={files} />}</Flex>
