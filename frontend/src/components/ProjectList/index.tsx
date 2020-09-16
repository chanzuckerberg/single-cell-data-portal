import React, { FC } from "react";
import { Project as IProject } from "src/common/entities";
import Dataset from "./components/Dataset";
import Heading from "./components/Heading";

interface Props {
  projects: IProject[];
}

const ProjectsList: FC<Props> = ({ projects }) => {
  if (!projects) return <div>Sorry, we could not find any projects</div>;

  return (
    <>
      <h1>Datasets</h1>
      <p>
        The cellxgene data portal is a repository of public, explorable
        single-cell datasets. If you have a public dataset which you would like
        hosted for visualization on this site, please let us know at{" "}
        <a href="mailto:cellxgene@chanzuckerberg.com">
          cellxgene@chanzuckerberg.com
        </a>
        .
      </p>
      <Heading />
      {projects
        .flatMap((project) =>
          project.datasets.map((dataset) => ({ dataset, links: project.links }))
        )
        .map(({ dataset, links }) => (
          <Dataset key={dataset.id} dataset={dataset} links={links} />
        ))}
    </>
  );
};

export default ProjectsList;
