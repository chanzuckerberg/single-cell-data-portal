import React, { FC } from "react";
import { useAuth0 } from "../contexts/auth0Context";
import { Text } from "theme-ui";

const LoginSignup: FC = () => {
  const auth = useAuth0();

  if (!auth) return null;

  const { isAuthenticated, loading, loginWithRedirect } = auth;

  if (loading || isAuthenticated) return null;

  return (
    <Text sx={{ color: "primary", fontSize: 1 }}>
      <a href="/" onClick={loginWithRedirect}>
        Log in
      </a>{" "}
      or{" "}
      <a href="/" onClick={loginWithRedirect}>
        sign-up
      </a>{" "}
      to view and download data.
    </Text>
  );
};

export default LoginSignup;
