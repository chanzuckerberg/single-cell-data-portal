import { Link } from "gatsby"
import React from "react"

const ProjectsList = ({ projects }) => {
  return (
    <div>
      {projects.map(project => {
        return (
          <div key={project.id}>
            <h3> {project.short_name} </h3>
            <Link to={`/project?id=${project.id}`}>{project.title}</Link>
            <p> {project.total_cells} cells </p>
            <p>
              Organs:
              {project.organs.map(organ => (
                <span key={organ}>{organ}</span>
              ))}
            </p>
            <p>
              Assays:
              {project.organs.map(assay => (
                <span key={assay}>{assay}</span>
              ))}
            </p>
            <p>
              Species:
              {project.organs.map(species => (
                <span key={species}>{species}</span>
              ))}
            </p>
          </div>
        )
      })}
    </div>
  )
}

export default ProjectsList

// short_name: "MouseGastrulationAtlas"
// title: "A single-cell molecular map of mouse gastrulation and early organogenesis"
// total_file_size: 339089
// id: "19"
// assays: ["10X 3' v1 sequencing"]
// total_cells: 339089
// organs: ["embryo"]
// species: ["Mus musculus"]
