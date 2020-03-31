import React from "react"
import { useAuth0 } from "../contexts/auth0Context"
import { Text } from "theme-ui"

const LoginSignup = () => {
  const { isAuthenticated, loading, loginWithRedirect } = useAuth0()

  return (
    <>
      {!loading && !isAuthenticated && (
        <Text sx={{ color: "primary", fontSize: 1 }}>
          <a href="/" onClick={loginWithRedirect}>Log in</a> or <a href="/" onClick={loginWithRedirect}>sign-up</a> to
          view and download data.
        </Text>)}
    </>
  )
}

export default LoginSignup
