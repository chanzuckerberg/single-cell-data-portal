import React, { FC } from "react";
import { useAuth0 } from "../contexts/auth0Context";
import { Box, Flex, Image } from "theme-ui";
import StyledButton from "./styledButton";
import { User } from "../common/entities";

interface LoggedInProps {
  user: User | null;
  logout: () => void;
}

const LoggedIn: FC<LoggedInProps> = ({ user, logout }) => {
  if (!user) return null;

  return (
    <Flex sx={{ alignItems: "center" }}>
      <StyledButton label="Logout" handleClick={logout} />
      <Image src={user.picture} width="45px" height="45px" padding="5px" />)
    </Flex>
  );
};

const Authenticate = () => {
  const auth = useAuth0();

  if (!auth) return null;

  const { isAuthenticated, loading, user, logout, loginWithRedirect } = auth;

  if (loading) return null;

  return (
    <Box>
      {<LoggedIn user={user} logout={logout} />}
      {!isAuthenticated && (
        <StyledButton label="Login/Sign-up" handleClick={loginWithRedirect} />
      )}
    </Box>
  );
};

export default Authenticate;
