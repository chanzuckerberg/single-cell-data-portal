import React, { FC, useEffect, useState } from "react";
import { Heading } from "theme-ui";
import { Project } from "../common/entities";
import Layout from "../components/layout";
import ProjectsList from "../components/projectsList";
import SEO from "../components/seo";
import { API_URL } from "../globals";

/*
  Mock API
  Get projects: https://ye54tu6ueg.execute-api.us-east-1.amazonaws.com/dev/projects
  Get project info: https://ye54tu6ueg.execute-api.us-east-1.amazonaws.com/dev/project/{id}
  Get project file: https://ye54tu6ueg.execute-api.us-east-1.amazonaws.com/dev/project/{id}/{file_name}
*/

const IndexPage: FC = () => {
  // Client-side Runtime Data Fetching
  const [projects, setProjects] = useState<Project[] | null>(null);

  useEffect(() => {
    fetch(`${API_URL}/projects`)
      .then(response => response.json()) // parse JSON from request
      .then(resultData => {
        setProjects(resultData);
      });
  }, []);

  return (
    <Layout>
      <SEO title="Explore Data" />
      <Heading
        as="h1"
        sx={{
          mb: 6,
          mt: 6,
        }}
      >
        Explore Data
      </Heading>
      {projects ? <ProjectsList projects={projects} /> : "Loading projects..."}
    </Layout>
  );
};

export default IndexPage;
