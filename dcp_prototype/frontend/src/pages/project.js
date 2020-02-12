import React, { useState, useEffect } from "react"

import Layout from "../components/layout"
import SEO from "../components/seo"
import searchStringAsObj from "../util/searchStringAsObj"
import FilesList from "../components/filesList"

const SecondPage = props => {
  const [project, setProject] = useState(null)
  const searchObj = searchStringAsObj(props.location.search.slice(1))
  const id = searchObj?.id

  useEffect(() => {
    if (id) {
      fetch(
        `https://ye54tu6ueg.execute-api.us-east-1.amazonaws.com/dev/project/${id}`
      )
        .then(response => response.json()) // parse JSON from request
        .then(resultData => {
          setProject(resultData)
        })
    } else {
      return null
    }
  }, [id]) /* don't rerender if this hasn't changed */

  return (
    <Layout>
      <SEO title="Projects" />
      <h1>Project</h1>
      <div>
        {!project?.files ? (
          "Loading project..."
        ) : (
          <FilesList files={project.files} />
        )}
      </div>
    </Layout>
  )
}

export default SecondPage
