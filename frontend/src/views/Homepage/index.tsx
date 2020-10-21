import React, { FC, useEffect, useState } from "react";
import { Project } from "src/common/entities";
import ProjectsList from "src/components/ProjectList";
import { API_URL } from "src/configs/configs";

const Homepage: FC = () => {
  const [projects, setProjects] = useState<Project[] | null>(null);
  const [isLoading, setIsLoading] = useState<boolean>(true);

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

export default Homepage;

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
    const response = await (await fetch(`${API_URL}/dp/v1/project`)).json();
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
  return (await fetch(`${API_URL}/dp/v1/project/${id}`)).json();
}
