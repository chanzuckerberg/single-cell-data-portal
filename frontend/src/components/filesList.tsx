import React, { FC } from "react";
import { Text, Flex, Box } from "theme-ui";
import { File } from "../common/entities";

interface Props {
  files: [File];
}

const FilesList: FC<Props> = ({ files }) => {
  console.log(files);
  return (
    <Box>
      {files.map(file => {
        return (
          <Flex key={file.id} sx={{ justifyContent: "space-between" }}>
            <Text sx={{ fontSize: 1 }}>Filename: {file.filename}</Text>
            <Text sx={{ fontSize: 1 }}>{file.species}</Text>
            <Text sx={{ fontSize: 1 }}>{file.file_size}</Text>
            <Text sx={{ fontSize: 1 }}>
              {file.library_construction_method_ontology}
            </Text>
            <Text sx={{ fontSize: 1 }}>{file.tissue_ontology}</Text>
            <Text sx={{ fontSize: 1 }}>{file.file_format}</Text>
          </Flex>
        );
      })}
    </Box>
  );
};

export default FilesList;
