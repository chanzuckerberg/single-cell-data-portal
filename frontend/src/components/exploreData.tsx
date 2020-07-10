import React, { FC } from "react";
import { Box, Flex, Text } from "theme-ui";
import { File, Project } from "../common/entities";
import { useAuth0 } from "../contexts/auth0Context";
import { API_URL, CXG_URL } from "../globals";
import LoginSignup from "./login-signup";
import StyledButton from "./styledButton";

interface Props {
  project: Project;
  files: [File];
  isAuthenticated: boolean;
}

const ExploreData: FC<Props> = ({ project, files, isAuthenticated }) => {
  const auth = useAuth0();

  const getTokenSilently = auth?.getTokenSilently;

  const matrixAvailable = (files && files.length && isAuthenticated) === true;
  const cxgAvailable = project.cxg_enabled;

  const handleClick = () => {
    if (!getTokenSilently) return;

    getTokenSilently()
      .then(accessToken =>
        fetch(`${API_URL}/files/${files[0].id}`, {
          headers: {
            Authorization: "Bearer " + accessToken,
          },
        })
      )
      .then(response => response.json())
      .then(resultData => {
        window.open(resultData.url);
      });
  };

  return (
    <Box
      sx={{
        mb: 2,
        maxWidth: 0,
      }}
    >
      <Flex
        sx={{
          mb: 1,
          maxWidth: 0,
        }}
      >
        <StyledButton
          label="Download matrix"
          disabled={!matrixAvailable}
          handleClick={handleClick}
        >
          Download matrix
        </StyledButton>
        <StyledButton
          label="Visualize in cellxgene"
          disabled={!cxgAvailable}
          handleClick={() => window.open(`${CXG_URL}/${project.label}.cxg`)}
        >
          Visualize in cellxgene
        </StyledButton>
      </Flex>
      {isAuthenticated ? null : (
        <Box>
          <LoginSignup />
        </Box>
      )}
      {files && !files.length && isAuthenticated ? (
        <Text sx={{ color: "primary", fontSize: 1 }}>
          No files available for download
        </Text>
      ) : null}
    </Box>
  );
};
export default ExploreData;
