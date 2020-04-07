import React, { useEffect, useState } from "react"
import { useAuth0 } from "../contexts/auth0Context"
import Layout from "../components/layout"
import SEO from "../components/seo"
import searchStringAsObj from "../util/searchStringAsObj"
import ProjectOverview from "../components/projectOverview"
import { api_url } from "../globals"
import { Flex, Heading } from "theme-ui"

const SecondPage = props => {
  const [project, setProject] = useState(null)
  const [files, setFiles] = useState(null)
  const searchObj = searchStringAsObj(props.location.search.slice(1))
  const id = searchObj?.id
  const { loading, isAuthenticated, getTokenSilently } = useAuth0()

  useEffect(() => {
    if (id) {
      fetch(`${api_url}/projects/${id}`)
        .then(response => response.json()) // parse JSON from request
        .then(resultData => {
          setProject(resultData)
        })
    }
  }, [id]) /* don't rerender if this hasn't changed */

  useEffect(() => {
    if (id && !loading && isAuthenticated) {
      // TODO send authorization token in requests to backend
      getTokenSilently().then(accessToken =>
        fetch(`${api_url}/projects/${id}/files`)
          .then(response => response.json()) // parse JSON from request
          .then(resultData => {
            setFiles(resultData)
          }))
    }
  }, [id, isAuthenticated, loading]) /* don't rerender if this hasn't changed */

  /*
    /projects/id
    /projects/id/files
    /files/id
  */

  return (
    <Layout>
      <SEO title="Projects"/>
      <Heading as="h1" sx={{ mb: 5, mt: 6 }}>
        Explore Project
      </Heading>
      <Flex sx={{justifyContent: 'center'}}>
        {!project ? (
          "Loading project..."
        ) : (
          <ProjectOverview project={project} files={files} isAuthenticated={isAuthenticated}/>
        )}
      </Flex>
    </Layout>
  )
}

export default SecondPage

// <Flex>{!files ? null : <FilesList files={files} />}</Flex>
