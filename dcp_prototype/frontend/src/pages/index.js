import React, { useState, useEffect } from "react"

import Layout from "../components/layout"
import SEO from "../components/seo"
import ProjectsList from "../components/projectsList"

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
    fetch(`https://ye54tu6ueg.execute-api.us-east-1.amazonaws.com/dev/projects`)
      .then(response => response.json()) // parse JSON from request
      .then(resultData => {
        setProjects(resultData)
      })
  }, [])
  return (
    <Layout>
      <SEO title="Home" />
      <h1>Explore Data</h1>
      {!projects ? "Loading projects..." : <ProjectsList projects={projects} />}
    </Layout>
  )
}

export default IndexPage
