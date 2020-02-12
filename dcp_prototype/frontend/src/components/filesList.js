import React from "react"

const FilesList = ({ files }) => {
  return (
    <div>
      {files.map(file => {
        return (
          <div key={file.file_key}>
            <h3>{file.file_key}</h3>
            <p>{file.species}</p>
            <p>{file.assay}</p>
            <p>{file.organ}</p>
            <p>{file.file_type}</p>
            <p>{file.file_name}</p>
            <p>{file.file_count}</p>
            <p>{file.file_size}</p>
            <a
              href={`https://ye54tu6ueg.execute-api.us-east-1.amazonaws.com/dev/project/${file.project_id}/${file.file_name}`}
              download={`${file.file_key}_${file.file_name}`}
            >
              Download File
            </a>
          </div>
        )
      })}
    </div>
  )
}

export default FilesList
