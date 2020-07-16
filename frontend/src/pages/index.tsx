import React, { FC, useEffect, useState } from "react";
import { API_URL } from "src/globals";
import { Project } from "../common/entities";
import Layout from "../components/Layout";
import ProjectsList from "../components/ProjectList";
import SEO from "../components/seo";

const Index: FC = () => {
  const [projects, setProjects] = useState<Project[] | null>(null);

  useEffect(() => {
    fetchProjects(setProjects);
  }, []);

  return (
    <Layout>
      <SEO title="Explore Data" />
      {projects ? <ProjectsList projects={projects} /> : "Loading projects..."}
    </Layout>
  );
};

interface ProjectResponse {
  project_uuid: string;
}

async function fetchProjects(setProjects: (allProjects: Project[]) => void) {
  const response = await (await fetch(`${API_URL}/v1/project`)).json();
  const projectIds: string[] = response.projects.map(
    (project: ProjectResponse) => project.project_uuid
  );

  const results = await Promise.allSettled(projectIds.map(fetchProject));

  const allProjects = results.reduce((memo, result) => {
    if (result.status === "fulfilled") {
      memo.push(result.value);
    }

    return memo;
  }, [] as Project[]);

  setProjects(allProjects);
}

async function fetchProject(id: string): Promise<Project> {
  return (await fetch(`${API_URL}/v1/project/${id}`)).json();
}

export default Index;
