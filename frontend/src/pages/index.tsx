import React, { FC, useEffect, useState } from "react";
import CookieBanner from "src/components/CookieBanner";
import ProjectsList from "src/components/ProjectList";
import { API_URL } from "src/configs/configs";
import { Project } from "../common/entities";
import Layout from "../components/Layout";
import SEO from "../components/seo";

const Content: FC = () => {
  const [projects, setProjects] = useState<Project[] | null>(null);
  const [isLoading, setIsLoading] = useState<boolean>(false);

  useEffect(() => {
    fetchProjects({
      setIsLoading,
      setProjects,
    });
  }, []);

  if (isLoading) return <div>Loading projects...</div>;

  return projects ? (
    <ProjectsList projects={projects} />
  ) : (
    <div>Could not find any projects. Please try again!</div>
  );
};

const Index: FC = () => {
  return (
    <Layout>
      <SEO title="Explore Data" />
      <Content />
      <CookieBanner />
    </Layout>
  );
};

interface ProjectResponse {
  id: string;
}

interface FetchProjectsArgs {
  setProjects: (allProjects: Project[]) => void;
  setIsLoading: (state: boolean) => void;
}

async function fetchProjects({ setProjects, setIsLoading }: FetchProjectsArgs) {
  setIsLoading(true);

  try {
    const response = await (await fetch(`${API_URL}/v1/project`)).json();
    const projectIds: string[] = response.projects.map(
      (project: ProjectResponse) => project.id
    );

    const results = await Promise.allSettled(projectIds.map(fetchProject));

    const allProjects = results.reduce((memo, result) => {
      if (result.status === "fulfilled") {
        memo.push(result.value);
      }

      return memo;
    }, [] as Project[]);

    setProjects(allProjects);
  } catch (error) {
    console.error(error);
  }

  setIsLoading(false);
}

async function fetchProject(id: string): Promise<Project> {
  return (await fetch(`${API_URL}/v1/project/${id}`)).json();
}

export default Index;
