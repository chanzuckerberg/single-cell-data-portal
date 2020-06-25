import React from "react";
import { useAuth0 } from "../contexts/auth0Context";
import { Box, Flex, Image } from "theme-ui";
import StyledButton from "./styledButton";

const Authenticate = () => {
  const auth = useAuth0();

  if (!auth) return null;

  const { isAuthenticated, loading, user, logout, loginWithRedirect } = auth;

  return (
    <Box>
      {!loading && user && (
        <Flex sx={{ alignItems: "center" }}>
          <StyledButton label="Logout" handleClick={logout} />
          {!user ? (
            <Box>...</Box>
          ) : (
            <Image
              src={user.picture}
              width="45px"
              height="45px"
              padding="5px"
            />
          )}
        </Flex>
      )}
      {!loading && !isAuthenticated && (
        <>
          <StyledButton label="Login/Sign-up" handleClick={loginWithRedirect} />
        </>
      )}
    </Box>
  );
};

export default Authenticate;
