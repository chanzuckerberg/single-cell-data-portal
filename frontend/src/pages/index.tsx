import React, { FC, useEffect, useState } from "react";
import { Project } from "../common/entities";
import { PROJECTS } from "../common/fixtures/projects";
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

async function fetchProjects(setProjects: (allProjects: Project[]) => void) {
  const projectIds: number[] = await new Promise(resolve => {
    setTimeout(() => {
      resolve(Array.from(Array(6)).map((_, index) => index));
    }, 0.5 * 1000);
  });

  const results = await Promise.allSettled(projectIds.map(fetchProject));

  const allProjects = results.reduce((memo, result) => {
    if (result.status === "fulfilled") {
      memo.push(result.value);
    }

    return memo;
  }, [] as Project[]);

  setProjects(allProjects);
}

async function fetchProject(id: number): Promise<Project> {
  return new Promise(resolve => {
    setTimeout(() => {
      resolve(PROJECTS[id]);
    }, 0.3 * 1000);
  });
}

export default Index;
