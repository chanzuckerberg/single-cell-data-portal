import React, { useState, useEffect } from "react"

import Layout from "../components/layout"
import SEO from "../components/seo"
import ProjectsList from "../components/projectsList"
import { Heading } from "theme-ui"
import { api_prefix } from "../globals"

/*
  Mock API
  Get projects: https://ye54tu6ueg.execute-api.us-east-1.amazonaws.com/dev/projects
  Get project info: https://ye54tu6ueg.execute-api.us-east-1.amazonaws.com/dev/project/{id}
  Get project file: https://ye54tu6ueg.execute-api.us-east-1.amazonaws.com/dev/project/{id}/{file_name}
*/

const IndexPage = () => {
  // Client-side Runtime Data Fetching
  const [projects, setProjects] = useState(null)
  useEffect(() => {
    fetch(`${api_prefix}/projects`)
      .then(response => response.json()) // parse JSON from request
      .then(resultData => {
        setProjects(resultData)
      })
  }, [])
  return (
    <Layout>
      <SEO title="Home" />
      <Heading as="h1" sx={{ mb: 4 }}>
        Explore Data
      </Heading>
      {!projects ? "Loading projects..." : <ProjectsList projects={projects} />}
    </Layout>
  )
}

export default IndexPage
